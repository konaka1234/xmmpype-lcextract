import os
import shutil
import numpy as np
from astropy.io import fits
from astropy.table import Table
import pxsas
import time

def extract_lc(obs_id, lc_bin=1000):
    # Set up directories
    work_dir = f"/data3/konakal/data/proc/{obs_id}/{obs_id}/"
    output_dir = f"/data3/konakal/data/lc/{obs_id}/"
    os.makedirs(output_dir, exist_ok=True)

    # Set the SAS_CCF environment variable to the CCF file in the obsid directory
    ccf_file = os.path.join(work_dir, 'ccf.cif')
    if os.path.exists(ccf_file):
        os.environ['SAS_CCF'] = ccf_file
    else:
        print(f"CCF file not found for OBSID {obs_id}. Proceeding without setting SAS_CCF.")

    # Find event file in the obsid directory
    eventfile = None
    for file_name in os.listdir(work_dir):
        if file_name.endswith("PIEVLI0000.FILTER"):
            eventfile = os.path.join(work_dir, file_name)
            break

    if eventfile is None:
        print(f"No event file found for OBSID {obs_id}.")
        return None

    # Iterate over all region files in the regions directory
    regions_dir = os.path.join(work_dir, 'regions')
    region_files = [f for f in os.listdir(regions_dir) if f.startswith(('src_', 'bkg_'))]

    # pair source and background regions
    region_dict = {}
    for region_file in region_files:
        if region_file.startswith('src_'):
            region_type = 'source'
        elif region_file.startswith('bkg_'):
            region_type = 'bkg'
        else:
            continue

        # Parse SDSS name from filename
        sdss_name = '_'.join(region_file.split('_')[1:-1])
        region_dict.setdefault(sdss_name, {})[region_type] = region_file

    for sdss_name, region_files in region_dict.items():
        source_lc_file = None
        bkg_lc_file = None

        # Extract source light curve if available
        if 'source' in region_files:
            region_file = region_files['source']
            region_path = os.path.join(regions_dir, region_file)

            # Read region coordinates from the region file
            with open(region_path, 'r') as f:
                region_data = f.read()

            # Extract coordinates from the source region file (x,y,radius)
            coordinates = region_data.split('(')[1].split(')')[0].split(',')
            x, y, r = coordinates[0], coordinates[1], coordinates[2]

            output_lc_file = f'{output_dir}{obs_id}_{sdss_name}_source.LC'

            # Define the filtering expression
            q_flag = "#XMMEA_EP"
            n_pattern = 4
            pn_pi_min = 500
            pn_pi_max = 2000
            expression = f"{q_flag}&&(PATTERN<={n_pattern})&&((X,Y) IN circle({x},{y},{r}))&&(PI in [{pn_pi_min}:{pn_pi_max}])"

            # Execute evselect to generate lc
            try:
                pxsas.run(
                    "evselect",
                    table=eventfile,
                    energycolumn="PI",
                    withrateset="yes",
                    rateset=output_lc_file,
                    timebinsize=lc_bin,
                    maketimecolumn="yes",
                    makeratecolumn="no",  # Ensure COUNTS column is used instead of RATE column
                    expression=expression
                )
                print(f"Generated light curve for OBSID {obs_id}, SDSS {sdss_name}, type: source")
                source_lc_file = output_lc_file
            except Exception as e:
                print(f"Failed to generate source light curve for OBSID {obs_id}, SDSS {sdss_name}. Error: {e}")

        # Extract background light curve if available
        if 'bkg' in region_files:
            region_file = region_files['bkg']
            region_path = os.path.join(regions_dir, region_file)

            # Read region coordinates from the region file
            with open(region_path, 'r') as f:
                region_data = f.read()

            # Extract coordinates from the background region file (x,y,inner,outer)
            coordinates = region_data.split('(')[1].split(')')[0].split(',')
            x, y, r_inner, r_outer = coordinates[0], coordinates[1], coordinates[2], coordinates[3]
            r = r_outer
            mask_file = os.path.join(work_dir, 'masks', region_file.replace('.reg', '.SRCMSK'))

            # Create temporary directory for the mask file
            temp_dir = os.path.join(output_dir, 'temp_mask')
            os.makedirs(temp_dir, exist_ok=True)
            temp_mask_file = os.path.join(temp_dir, 'bkg.SRCMSK')

            # Copy mask file to temporary directory with simple name
            try:
                shutil.copy(mask_file, temp_mask_file)
            except Exception as e:
                print(f"Failed to copy mask file {mask_file} to temporary directory. Error: {e}")
                continue

            output_lc_file = f'{output_dir}{obs_id}_{sdss_name}_bkg.LC'

            # Define filtering expression
            q_flag = "#XMMEA_EP"
            n_pattern = 4
            pn_pi_min = 500
            pn_pi_max = 2000
            expression = (
                f"{q_flag}&&(PATTERN<={n_pattern})"
                f"&&mask({temp_mask_file},0,0,X,Y)&&(PI in [{pn_pi_min}:{pn_pi_max}])&&(X,Y) in annulus({x},{y},{r_inner},{r})"
            )

            # Execute evselect
            try:
                pxsas.run(
                    "evselect",
                    table=eventfile,
                    energycolumn="PI",
                    withrateset="yes",
                    rateset=output_lc_file,
                    timebinsize=lc_bin,
                    maketimecolumn="yes",
                    makeratecolumn="no",  # Ensure COUNTS column is used instead of RATE column
                    expression=expression
                )
                print(f"Generated light curve for OBSID {obs_id}, SDSS {sdss_name}, type: bkg")
                bkg_lc_file = output_lc_file
            except Exception as e:
                print(f"Failed to generate background light curve for OBSID {obs_id}, SDSS {sdss_name}. Error: {e}")
            finally:
                # Clean up the temp directory
                try:
                    shutil.rmtree(temp_dir)
                except Exception as e:
                    print(f"Failed to remove temporary directory {temp_dir}. Error: {e}")

        # Generate corrected light curve if both source and background light curves are available
        if source_lc_file and bkg_lc_file:
            # Create a temp directory for the light curve files
            temp_lc_dir = os.path.join(output_dir, 'temp_lc')
            os.makedirs(temp_lc_dir, exist_ok=True)

            # Copy source and background light curves to the temp directory with simpler names
            temp_source_lc_file = os.path.join(temp_lc_dir, 'source.LC')
            temp_bkg_lc_file = os.path.join(temp_lc_dir, 'bkg.LC')
            try:
                shutil.copy(source_lc_file, temp_source_lc_file)
                shutil.copy(bkg_lc_file, temp_bkg_lc_file)
            except Exception as e:
                print(f"Failed to copy light curve files to temporary directory for corrected light curve. Error: {e}")
                continue

            corrected_lc_file = f'{output_dir}{obs_id}_{sdss_name}_corrlc.LC'
            try:
                pxsas.run(
                    "epiclccorr",
                    srctslist=temp_source_lc_file,
                    eventlist=eventfile,
                    outset=corrected_lc_file,
                    bkgtslist=temp_bkg_lc_file,
                    withbkgset="yes",
                    applyabsolutecorrections="yes"
                )
                print(f"Generated corrected light curve for OBSID {obs_id}, SDSS {sdss_name}")

                # Remove rows with FRACEXP v of 0 or NULL from corrected light curve
                with fits.open(corrected_lc_file, mode='update') as hdul:
                    lc_data = Table(hdul[1].data)
                    valid_rows = lc_data['FRACEXP'] > 0
                    filtered_data = lc_data[valid_rows]
                    hdul[1].data = filtered_data.as_array()
                    print(f"Filtered out rows with FRACEXP = 0 or NULL in corrected light curve for OBSID {obs_id}, SDSS {sdss_name}")
            except Exception as e:
                print(f"Failed to generate corrected light curve for OBSID {obs_id}, SDSS {sdss_name}. Error: {e}")
            finally:
                # Clean up the temporary directory after processing the corrected light curve
                try:
                    shutil.rmtree(temp_lc_dir)
                except Exception as e:
                    print(f"Failed to remove temporary directory {temp_lc_dir}. Error: {e}")

if __name__ == "__main__":
    test_obsid = "0693540401" 
    extract_lc(test_obsid)

