import logging
import csv
import xmmpype as xmm
from xmmpype.crossmatch import XMatch
import xmmpype.hpixels as xmmhp
from xmmpype.obsids import XMMPYobsid
from makesrclist import process_healpix_cells
from makereg import make_ds9regions
from makeqsoreg import generate_qso_regions
from excludesources import exclude_regions, create_sources_mask
from makebkgmask import create_bkg_masks
from corrlc import extract_lc
import os
import shutil
from multiprocessing import Pool
import time

# Function to process each obsid independently
def process_obsid(obsid):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    start_time = time.time()  # Start timing
    try:
        # Define project
        project_name = f"{obsid}"
        P = xmm.Project(root_folder='/home/konaka/xmmpype_extend/data/', mergedir='hp',
                        proc='proc', raw='raw', project_name=project_name,
                        astrocor_survey=None, eband="all", dbfile=f"{project_name}.db")

        # Add the single obsid to the project
        P.add_obsids([obsid])
        
        # Reduce the obsid
        P.reduce_obsids(ncores=2)

        # Define the HEALPix grid for the single obsid
        hpixels = P.calc_hpixels()

        # Print the identification numbers of the cells of the grid
        logger.info(f"HEALPix cells for OBS_ID {obsid}: {P.hpixels}")

        # Reduce hpixels (merge if there are overlapping observations, won't happen since we are using 1)
        P.reduce_hpixels(ncores=2)

        # Define the MOC of the Project
        P.mastermoc()

        # Define a list of unique sources for the project
        P.srclist()

        # Calculate the sensitivity map for the project
        P.sensemap()

        # Extract counts (make source lists) from healpix cells
        process_healpix_cells(project_name=obsid, eef=70)

        # Make ds9 regions for all sources in obsid
        make_ds9regions(obsid)

        # Make source and background regions for QSOs (from my catalog)
        generate_qso_regions(obsid)

        # Make a file masking out all sources
        create_sources_mask(obsid)

        # Make masks that mask out everything but background annulus for each source (also all other sources from the above)
        create_bkg_masks(obsid)

        # Extract Source, Background and Corrected lightcurves for each source
        extract_lc(obsid, lc_bin=1000)

        # Log successful processing
        logger.info(f"Completed processing for OBS_ID: {obsid}")

        # Move .log and .db files to the appropriate directories
        destination_log_directory = '/home/konaka/xmmpype_extend/logs/'
        destination_db_directory = '/home/konaka/xmmpype_extend/db/'
        os.makedirs(destination_log_directory, exist_ok=True)
        os.makedirs(destination_db_directory, exist_ok=True)

        log_files = [f for f in os.listdir('.') if f.endswith('.log')]
        for log_file in log_files:
            shutil.move(log_file, os.path.join(destination_log_directory, log_file))

        db_file = f"{project_name}.db"
        if os.path.exists(db_file):
            shutil.move(db_file, os.path.join(destination_db_directory, db_file))

        # Delete unnecessary raw files to save space
        raw_dir = f"/home/konaka/xmmpype_extend/data/raw/{obsid}/{obsid}/ODF/"
        if os.path.exists(raw_dir):
            for filename in os.listdir(raw_dir):
                if not (filename.endswith("SUM.ASC") or filename.endswith("SUM.SAS")):
                    os.remove(os.path.join(raw_dir, filename))

        end_time = time.time()  # End timing
        elapsed_time = end_time - start_time
        logger.info(f"Processing time for OBS_ID {obsid}: {elapsed_time:.2f} seconds")

        return obsid, True, elapsed_time  # Return success and elapsed time

    except Exception as e:
        logger.error(f"Failed processing for OBS_ID: {obsid} with error: {e}")
        return obsid, False, None  # Return failure

# Function to get obsids from QSO CSV 
def get_obsids_from_csv(file_path, max_obsids=None):
    obsids = set()
    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            obsids.add(row['OBS_ID'])
            if max_obsids is not None and len(obsids) >= max_obsids:
                break
    return [str(obsid) for obsid in obsids]  # Convert to string format

if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(level=logging.INFO, filename='/home/konaka/xmmpype_extend/logs/pipeline.log', filemode='w',
                        format='%(asctime)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger('').addHandler(console)

    # Load obsids from CSV file and add to the project
    csv_file = '/home/konaka/xmmpype_extend/data/catalogs/qso_coords.csv'
    max_obsids = 5 # Set to None to grab all unique values

    # Manual obsids to be added
    manual_obsids = [""]

    if max_obsids == 0:
        obsids_from_csv = manual_obsids
    else:
        obsids_from_csv = get_obsids_from_csv(csv_file, max_obsids) + manual_obsids 

    # Remove duplicates that exist in both lists
    obsids_from_csv = list(set(obsids_from_csv))

    start_time_total = time.time()  # Start total processing time

    try:
        pool = Pool(processes=4)
        results = pool.map(process_obsid, obsids_from_csv)
    except KeyboardInterrupt:
        # Close all running processes if CTRL+C is pressed
        pool.close()
        pool.join()
    finally:
        # To make sure processes are closed in the end, even if errors happen
        pool.close()
        pool.join()

    # Tracking processed, failed OBS_IDs and calculating average time
    processed_obsids = [obsid for obsid, success, _ in results if success]
    failed_obsids = [obsid for obsid, success, _ in results if not success]
    times = [time for _, success, time in results if success]
    avg_time = sum(times) / len(times) if times else 0

    total_obsids = len(obsids_from_csv)
    processed_count = len(processed_obsids)
    failed_count = len(failed_obsids)
    remaining_count = total_obsids - processed_count - failed_count

    end_time_total = time.time()  # End total processing time
    total_elapsed_time = end_time_total - start_time_total

    logging.info("\nProgress Update:")
    logging.info(f"Total OBS_IDs: {total_obsids}")
    logging.info(f"Successfully Processed: {processed_count}")
    logging.info(f"Failed: {failed_count}")
    logging.info(f"Remaining: {remaining_count}")

    # Final summary
    logging.info("\nProcessing Summary:")
    logging.info(f"Total OBS_IDs: {total_obsids}")
    logging.info(f"Successfully Processed: {processed_count}")
    logging.info(f"Failed: {failed_count}")
    logging.info(f"Remaining: {remaining_count}")
    logging.info(f"Average processing time per OBS_ID: {avg_time:.2f} seconds")
    logging.info(f"Total processing time for all OBS_IDs: {total_elapsed_time:.2f} seconds")

    if failed_count > 0:
        logging.info("\nFailed OBS_IDs:")
        for failed_obsid in failed_obsids:
            logging.info(f" - {failed_obsid}")

    logging.info("All OBS_IDs have been processed.")

