from astropy.table import Table
from astropy.io import fits
import os
from xmmpype.utils.coordinates import sky2phys  

def make_ds9regions(obsid):
    
    obsid_directory = f'/data3/konakal/data/proc/{obsid}/{obsid}'
    hp_directory = f'/data3/konakal/data/hp/{obsid}'
    
    # Path to the counts FITS file
    src_path = os.path.join(hp_directory, 'SRC', 'extracted_counts.fits')

    if os.path.exists(src_path):
        
        counts_table = Table.read(src_path)

        # Filter out rows with -100 or -9 values?
        valid_counts_table = counts_table[(counts_table['CNT'] != -100) &
                                          (counts_table['BKG'] != -100) &
                                          (counts_table['SRC_MAX'] != -9) &
                                          (counts_table['SRC_MEAN'] != -9)]

        ra_values = valid_counts_table['RA']
        dec_values = valid_counts_table['DEC']
        radii_arcsec = valid_counts_table['RADIUS']

        # Use the file ending in PIEVLI0000_FULL.IMG for WCS transformation
        img_file = [f for f in os.listdir(obsid_directory) if f.endswith("PIEVLI0000_FULL.IMG")][0]
        img_fits_path = os.path.join(obsid_directory, img_file)

        # Convert RA, DEC to detector coordinates (physical coordinates)
        x, y, rphys = sky2phys(ra_values, dec_values, wcs=None, img=img_fits_path, r=radii_arcsec)

        # Write the DS9 regions
        output_reg_file = os.path.join(obsid_directory, f'ds9_regions_{obsid}.reg')
        with open(output_reg_file, 'w') as reg_file:
            reg_file.write('physical\n')
            for xi, yi, radius in zip(x, y, rphys):
                reg_file.write(f'circle({xi},{yi},{radius})\n')

        print(f"DS9 region file created for OBSID {obsid}: {output_reg_file}")
    else:
        print(f"Counts FITS file not found for OBSID {obsid}")

if __name__ == "__main__":
    test_obsid = "0201900101"  
    make_ds9regions(test_obsid)

