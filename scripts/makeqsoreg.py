import os
import csv
from astropy.io import fits
from astropy.wcs import WCS
from xmmpype.utils.coordinates import sky2phys  

# Set the fixed scaling factor: 1 pixel = 4 arcseconds
arcsec_per_pixel = 4.0

def generate_qso_regions(obsid):
    
    obsid_directory = f'/data3/konakal/data/proc/{obsid}/{obsid}'
    hp_directory = f'/data3/konakal/data/hp/{obsid}'
    qso_catalog = '/data3/konakal/data/catalogs/qso_coords_new.csv'

    # Use the file ending in PIEVLI0000_FULL.IMG for WCS transformation
    try:
        img_file = [f for f in os.listdir(obsid_directory) if f.endswith("PIEVLI0000_FULL.IMG")][0]
    except IndexError:
        print(f"No appropriate image file found for OBSID: {obsid}")
        return

    img_fits_path = os.path.join(obsid_directory, img_file)

    # Read the QSO catalog csv
    with open(qso_catalog, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        qso_list = list(reader)

    success_count = 0
    not_found_count = 0

    # Process each QSO in the catalog
    for qso_data in qso_list:
        if qso_data['OBS_ID'] == obsid:
            ra = float(qso_data['RA'])
            dec = float(qso_data['DEC'])
            sdss_name = qso_data['SDSS_NAME']

            # Convert RA, DEC to physical detector coordinates using sky2phys
            x_pix, y_pix, _ = sky2phys([ra], [dec], wcs=None, img=img_fits_path, r=[arcsec_per_pixel])

            # Create the regions directory
            regions_dir = os.path.join(obsid_directory, 'regions')
            os.makedirs(regions_dir, exist_ok=True)

            # Read the main region file to get the radius for each source
            main_region_file = os.path.join(obsid_directory, f'ds9_regions_{obsid}.reg')
            radius_pix = None

            if os.path.exists(main_region_file):
                with open(main_region_file, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if line.startswith('circle'):
                            # Parse the region line to get the X, Y, and radius in physical coordinates
                            parts = line.strip().split('(')[1].split(')')[0].split(',')
                            x_region, y_region, radius_region = map(float, parts[:3])
                            # Compare the physical coordinates 
                            if abs(x_region - x_pix[0]) < 100 and abs(y_region - y_pix[0]) < 100:  # I should look into this more - used 100 pixels for now
                                radius_pix = radius_region
                                break
            else:
                print(f"No region file found for OBSID: {obsid}")
                not_found_count += 1
                continue

            if radius_pix is None:
                print(f"No matching region found for QSO at RA: {ra}, DEC: {dec} in OBSID: {obsid}")
                print(f"Physical coordinates: X = {x_pix[0]}, Y = {y_pix[0]}")
                not_found_count += 1
                continue

            # Define inner and outer radii for the background annulus in pixels
            inner_radius_pix = radius_pix
            outer_radius_pix = radius_pix + 2480

            # Define paths for source and background region files
            source_region_file = os.path.join(regions_dir, f'src_{sdss_name}_{obsid}.reg')
            background_region_file = os.path.join(regions_dir, f'bkg_{sdss_name}_{obsid}.reg')

            # Save the source region
            with open(source_region_file, 'w') as f:
                f.write('physical\n')
                f.write(f'circle({x_pix[0]},{y_pix[0]},{radius_pix}) # color=white\n')

            # Save the background region
            with open(background_region_file, 'w') as f:
                f.write('physical\n')
                f.write(f'annulus({x_pix[0]},{y_pix[0]},{inner_radius_pix},{outer_radius_pix}) # color=magenta\n')

            success_count += 1

    print(f"Finished processing OBSID {obsid}. Success: {success_count}, Not Found: {not_found_count}")

if __name__ == "__main__":
    test_obsid = "0201900101"  
    generate_qso_regions(test_obsid)

