import os
import numpy as np
from astropy.io import fits
from xmmpype.utils.coordinates import phys2pix 


def create_bkg_masks(obsid):
   
    obsid_dir = f'/data3/konakal/data/proc/{obsid}/{obsid}'
    base_image_path = os.path.join(obsid_dir, f'{obsid}.TOTALSRCMSK')
    filter_image_candidates = [f for f in os.listdir(obsid_dir) if f.endswith("PIEVLI0000_FULL.IMG")]
    if not filter_image_candidates:
        raise FileNotFoundError(f"No .IMG file found for WCS transformation in OBSID {obsid}")
    filter_image_path = os.path.join(obsid_dir, filter_image_candidates[0])

   
    with fits.open(base_image_path, mode='readonly', output_verify="silentfix") as hdu:
        original_maskdata = hdu[1].data

   
    masks_dir = os.path.join(obsid_dir, 'masks')
    os.makedirs(masks_dir, exist_ok=True)

   
    regions_dir = os.path.join(obsid_dir, 'regions')
    region_files = [f for f in os.listdir(regions_dir) if f.startswith('bkg')]

    for region_file in region_files:
       
        maskdata = np.zeros_like(original_maskdata)

        
        region_path = os.path.join(regions_dir, region_file)
        with open(region_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("annulus"):
                    
                    parts = line.split('(')[1].split(')')[0].split(',')
                    x_phys, y_phys, inner_radius_phys, outer_radius_phys = map(float, parts[:4])

                    # Convert physical coordinates to pixel coordinates
                    try:
                        x_pix, y_pix, inner_radius_pix = phys2pix(x_phys, y_phys, img=filter_image_path, rphys=inner_radius_phys)
                        _, _, outer_radius_pix = phys2pix(x_phys, y_phys, img=filter_image_path, rphys=outer_radius_phys)
                    except Exception as e:
                        print(f"Error converting physical coordinates to pixel coordinates: {e}")
                        continue

                    # Set only the pixels within the annulus to their original values
                    for y in range(maskdata.shape[0]):
                        for x in range(maskdata.shape[1]):
                            dist = np.sqrt((x - x_pix) ** 2 + (y - y_pix) ** 2)
                            if inner_radius_pix < dist <= outer_radius_pix:
                                maskdata[y, x] = original_maskdata[y, x]

        # Save the mask data for region file WITHOUT altering the original file
        output_mask_path = os.path.join(masks_dir, f"{region_file.replace('.reg', '.SRCMSK')}")
        hdu_mask = fits.PrimaryHDU(data=maskdata, header=hdu[1].header)
        hdu_mask.writeto(output_mask_path, overwrite=True)
        print(f"Mask created successfully: {output_mask_path}")


if __name__ == "__main__":
    test_obsid = "0201900101"  
    create_bkg_masks(test_obsid)

