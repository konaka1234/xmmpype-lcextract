import logging
from astropy.table import Table, vstack
import xmmpype as xmm
import xmmpype.hpixels as xmmhp
from xmmpype.obsids import XMMPYobsid
from xmmpype.events import XMMPYevt
import os

logging.basicConfig(level=logging.INFO)

def process_healpix_cells(root_folder='/data3/konakal/data/', project_name='Test1', eef=70):
    """
    Process HEALPix cells to extract counts and write to FITS files.

    Parameters:
    - root_folder (str): Root folder of the project data.
    - project_name (str): Name of the project.
    - eef (int): Effective extraction fraction.
    """
    # Initialize the project
    P = xmm.Project(root_folder=root_folder, mergedir='hp',
                    proc='proc', raw='raw', project_name=project_name,
                    astrocor_survey=None, eband="all", dbfile="{}.db".format(project_name))

    # Loop over all HEALPix directories (assuming they are named with numbers)
    healpix_dirs = [d for d in os.listdir(os.path.join(root_folder, 'hp', project_name)) if d.isdigit()]

    extracted_counts_tables = []

    for hpix_dir in healpix_dirs:
        src_path = os.path.join(root_folder, 'hp', project_name, hpix_dir, 'SRC', 'srclist.fits')

        if os.path.exists(src_path):
            logging.info("Processing srclist.fits in HEALPix cell {}...".format(hpix_dir))

            # Load the srclist.fits file
            srclist_table = Table.read(src_path, hdu=1)

            # Extract the RA and DEC columns
            srclist = Table()
            srclist['RA'] = srclist_table['RA']
            srclist['DEC'] = srclist_table['DEC']

            # Initialize HEALPix instance for the current HEALPix cell
            hp = xmmhp.HEALpix(P, int(hpix_dir))
            eband = hp.project.ebands[1]

            for evtf in hp.event_files:
                if evtf['detector'] == "PN":
                    obs = XMMPYobsid(hp.project, evtf["obsid"])
                    evt = XMMPYevt(obs, evtf["filename"])
                    counts = evt.extract_counts(eband, srclist, eef, bkg_path=hp.paths.tmp)
                    logging.info("Counts extracted for HEALPix cell {}: {}".format(hpix_dir, counts))

                    # Define output path for extracted counts
                    output_fits_file = os.path.join(root_folder, 'hp', project_name, hpix_dir, 'SRC', 'extracted_counts.fits')

                    # Overwrite the FITS file with the new data
                    counts.write(output_fits_file, format='fits', overwrite=True)
                    logging.info("Data written to FITS file: {}".format(output_fits_file))

                    # Add the extracted counts to the list for combining later
                    extracted_counts_tables.append(counts)

                    logging.info("HEALPIX Cell: {}; BAND: {}; OBSID: {}; FILENAME: {}; DETECTOR: {}".format(
                        hpix_dir, eband, evtf['obsid'], evtf['filename'], evtf['detector']))
        else:
            logging.info("No srclist.fits found in HEALPix cell {}, skipping.".format(hpix_dir))

    # Combine all extracted counts tables into one and write to the total SRC directory
    if extracted_counts_tables:
        combined_counts = vstack(extracted_counts_tables)
        combined_output_path = os.path.join(root_folder, 'hp', project_name, 'SRC', 'extracted_counts.fits')
        combined_counts.write(combined_output_path, format='fits', overwrite=True)
        logging.info("Combined extracted counts written to: {}".format(combined_output_path))

    logging.info("Processing complete for all HEALPix cells.")


def main():
    obsid = "0201900101"  
    process_healpix_cells(project_name=obsid)

if __name__ == "__main__":
    main()

