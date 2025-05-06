import os
import numpy as np
from astropy.io import fits
from xmmpype.utils.coordinates import phys2pix  
from regions import CirclePixelRegion
from regions.core import PixCoord
from astropy.wcs import WCS


def exclude_regions(image_data, regions):
    # Create a total mask to store all regions combined
    total_srcmask = np.zeros_like(image_data, dtype=float)
    
    for region in regions:
        # Ensure region is treated as a physical coordinate region
        if isinstance(region, CirclePixelRegion):
            mask = region.to_mask(mode='exact')  # Use 'exact' mode to cover the entire region hopefuly
            mask_data = mask.to_image(image_data.shape)
            if mask_data is not None:
                mask_data = np.nan_to_num(mask_data, nan=0.0)  # Replace NaNs with 0 to avoid issues (areas not in detector?)
                total_srcmask += mask_data  # Accumulate the mask for all regions
    
    # Set all pixels in the regions to 0
    image_data[total_srcmask > 0] = 0
    return image_data


def create_sources_mask(obsid):
    # Define paths 
    obsid_directory = f'/data3/konakal/data/proc/{obsid}/{obsid}'
    mask_file_path = os.path.join(obsid_directory, [f for f in os.listdir(obsid_directory) if f.endswith("PIEVLI0000_FULL.MSK")][0])
    regions_file_path = os.path.join(obsid_directory, f'ds9_regions_{obsid}.reg')
    filter_image_candidates = [f for f in os.listdir(obsid_directory) if f.endswith("PIEVLI0000_FULL.IMG")]
    if not filter_image_candidates:
        return
    filter_image_path = os.path.join(obsid_directory, filter_image_candidates[0])
    
    # Load the mask FITS file
    with fits.open(mask_file_path, mode='update', output_verify="silentfix") as hdu:
        maskhdr = hdu[1].header
        maskdata = hdu[1].data

        # Read the regions file to get regions of all sources
        if os.path.exists(regions_file_path):
            regions = []
            with open(regions_file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('circle'):
                        # Parse the region line to get the physical coordinates and radius
                        parts = line.split('(')[1].split(')')[0].split(',')
                        x_phys, y_phys, radius_phys = map(float, parts)
                        
                        # Convert physical detector coordinates to pixel coordinates using phys2pix
                        try:
                            x_pix, y_pix, rpix = phys2pix(x_phys, y_phys, wcs=None, img=filter_image_path, rphys=radius_phys)
                        except Exception as e:
                            continue

                        center = PixCoord(x_pix, y_pix)
                        region = CirclePixelRegion(center=center, radius=rpix)
                        regions.append(region)
            
            # Exclude the regions from the image through the function
            new_maskdata = exclude_regions(maskdata, regions)

            # Save the updated mask file with the sources masked out
            new_mask_file_path = os.path.join(obsid_directory, f'{obsid}.TOTALSRCMSK')
            hdu[1].data = new_maskdata
            hdu.writeto(new_mask_file_path, overwrite=True)
            print(f"Mask file with sources removed saved for OBSID {obsid}: {new_mask_file_path}")
        else:
            pass


if __name__ == "__main__":
    test_obsid = "0201900101"  
    create_sources_mask(test_obsid)

