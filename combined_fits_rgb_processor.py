#!/usr/bin/env python3
"""
Combined JWST FITS Preprocessing and RGB Image Creation Pipeline

This script combines FITS preprocessing and RGB image creation into a single workflow.
It processes entire directory structures, creating both preprocessed FITS files and RGB images.

Requirements:
    pip install astropy reproject matplotlib numpy Pillow

Usage:
    python combined_fits_rgb_processor.py --input-dir ./data --output-dir ./processed_fits_rgb
"""

import numpy as np
import argparse
import warnings
import os
import tempfile
import shutil
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from reproject import reproject_interp
import re

# Check for required Pillow library
try:
    from PIL import Image
    PILLOW_AVAILABLE = True
except ImportError:
    PILLOW_AVAILABLE = False
    print("WARNING: Pillow not found. Please install: pip install Pillow")

def apply_ds9_normalization(data, contrast=1.02, bias=0.43, max_percent=99):
    """Apply DS9-style normalization with asinh stretch and color inversion."""
    # Handle NaN and infinite values
    finite_mask = np.isfinite(data)
    if not np.any(finite_mask):
        print("Warning: No finite values in data")
        return data
    
    finite_data = data[finite_mask]
    
    # Calculate statistics on finite data
    mean, median, std = sigma_clipped_stats(finite_data, sigma=3.0)
    
    # Calculate min and max values
    min_val = np.min(finite_data)
    max_val = np.percentile(finite_data, max_percent)
    
    print(f"      Data stats: min={min_val:.3f}, max={max_val:.3f}, median={median:.3f}, std={std:.3f}")
    
    # Normalize to [0, 1] range first
    data_range = max_val - min_val
    if data_range == 0:
        print("Warning: Data has no dynamic range")
        return np.zeros_like(data)
    
    normalized = (data - min_val) / data_range
    normalized = np.clip(normalized, 0, 1)
    
    # Apply DS9 asinh stretch with contrast and bias
    stretched = np.arcsinh(contrast * (normalized - bias)) / np.arcsinh(contrast)
    
    # Invert colors (black becomes white, white becomes black)
    stretched = 1.0 - stretched
    
    # Restore NaN/inf values
    stretched[~finite_mask] = data[~finite_mask]
    
    return stretched

def preprocess_jwst_fits(input_path, output_path, contrast=1.02, bias=0.43, max_percent=99):
    """Convert JWST FITS file with SCI extension to single extension format with normalization."""
    try:
        # Open the FITS file
        with fits.open(input_path) as hdul:
            print(f"    Processing: {os.path.basename(input_path)}")
            
            # Find the SCI extension
            sci_ext = None
            for i, hdu in enumerate(hdul):
                if hdu.name == 'SCI' or (hasattr(hdu, 'header') and 
                                        hdu.header.get('EXTNAME', '').upper() == 'SCI'):
                    sci_ext = i
                    break
            
            if sci_ext is None:
                # If no SCI extension found, try the first extension with data
                for i, hdu in enumerate(hdul):
                    if hasattr(hdu, 'data') and hdu.data is not None:
                        sci_ext = i
                        break
            
            if sci_ext is None:
                raise ValueError("No suitable data extension found")
            
            # Extract header and data from SCI extension
            header = hdul[sci_ext].header.copy()
            data = hdul[sci_ext].data.copy()
            
            if data is None:
                raise ValueError(f"No data found in extension {sci_ext}")
            
            print(f"      Data shape: {data.shape}, dtype: {data.dtype}")
            
            # Apply DS9 normalization
            print(f"      Applying DS9 normalization with color inversion")
            normalized_data = apply_ds9_normalization(data, contrast, bias, max_percent)
            
            # Update header to indicate processing
            header['HISTORY'] = f'Processed with DS9 normalization: contrast={contrast}, bias={bias}'
            header['HISTORY'] = f'Applied asinh stretch with {max_percent}% max percentile'
            header['HISTORY'] = 'Applied color inversion (black <-> white)'
            header['DSCONTRST'] = (contrast, 'DS9 contrast parameter used')
            header['DSBIAS'] = (bias, 'DS9 bias parameter used')
            header['DSMAXPCT'] = (max_percent, 'Max percentile used for normalization')
            header['INVERTED'] = (True, 'Colors inverted (black <-> white)')
            
            # Remove extension-specific keywords that might cause issues
            remove_keywords = ['EXTNAME', 'EXTVER', 'INHERIT']
            for keyword in remove_keywords:
                if keyword in header:
                    del header[keyword]
            
            # Create new HDU with normalized data
            primary_hdu = fits.PrimaryHDU(data=normalized_data, header=header)
            
            # Create HDU list and write to file
            new_hdul = fits.HDUList([primary_hdu])
            
            # Create output directory if it doesn't exist
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            # Write the file
            new_hdul.writeto(output_path, overwrite=True)
            print(f"      Saved preprocessed FITS: {os.path.basename(output_path)}")
            
            return True
            
    except Exception as e:
        print(f"    Error processing {input_path}: {str(e)}")
        return False

def load_fits_with_wcs(filename):
    """Load FITS data and WCS information."""
    try:
        with fits.open(filename) as hdul:
            # Use primary HDU since we preprocessed to single extension
            data_hdu = hdul[0]
            
            if data_hdu.data is None:
                raise ValueError("No data found in primary HDU")
            
            data = data_hdu.data.copy()
            header = data_hdu.header.copy()
            
            try:
                wcs = WCS(header)
                # Handle pixel scale units properly
                try:
                    scales = wcs.proj_plane_pixel_scales()
                    if hasattr(scales[0], 'unit'):
                        if scales[0].unit.to_string() == 'deg':
                            scale_arcsec = scales[0].to_value('arcsec')
                        else:
                            scale_arcsec = scales[0].value * 3600
                    else:
                        scale_arcsec = float(scales[0]) * 3600
                except Exception:
                    scale_arcsec = 0.1  # fallback
            except Exception as e:
                print(f"      Warning: Could not load WCS: {e}")
                wcs = None
                scale_arcsec = 0.1
            
            return data, header, wcs, scale_arcsec
            
    except Exception as e:
        print(f"    Error loading {filename}: {e}")
        raise

def find_common_wcs(wcs_list, data_shapes, pixel_scales):
    """Determine the best common WCS grid for reprojection."""
    valid_wcs = [(wcs, shape, scale) for wcs, shape, scale in zip(wcs_list, data_shapes, pixel_scales) if wcs is not None]
    
    if not valid_wcs:
        # If no WCS available, use largest image as reference
        max_pixels = max(np.prod(shape) for shape in data_shapes)
        ref_idx = next(i for i, shape in enumerate(data_shapes) if np.prod(shape) == max_pixels)
        return None, data_shapes[ref_idx]
    
    # Use finest resolution as reference
    ref_idx = np.argmin([scale for _, _, scale in valid_wcs])
    ref_wcs = valid_wcs[ref_idx][0]
    ref_shape = valid_wcs[ref_idx][1]
    
    return ref_wcs, ref_shape

def reproject_to_common_grid(data_list, wcs_list, data_shapes, pixel_scales, ref_wcs, ref_shape):
    """Reproject all images to a common WCS grid."""
    reprojected_data = []
    
    for i, (data, wcs, shape, scale) in enumerate(zip(data_list, wcs_list, data_shapes, pixel_scales)):
        if ref_wcs is None or wcs is None:
            # No WCS available, use simple resizing
            if data.shape != ref_shape:
                reprojected = resize_to_shape(data, ref_shape)
            else:
                reprojected = data.copy()
        else:
            try:
                # Reproject to common WCS grid
                reprojected, footprint = reproject_interp(
                    (data, wcs), ref_wcs, shape_out=ref_shape,
                    order='bilinear'
                )
                
                # Handle NaN values from reprojection
                nan_mask = ~np.isfinite(reprojected)
                if np.any(nan_mask):
                    reprojected[nan_mask] = 0
                
            except Exception as e:
                print(f"      Warning: Reprojection failed ({e}), using resize fallback")
                reprojected = resize_to_shape(data, ref_shape)
        
        reprojected_data.append(reprojected)
    
    return reprojected_data

def resize_to_shape(data, target_shape):
    """Simple resize/crop to match target shape."""
    current_shape = data.shape
    target_h, target_w = target_shape
    current_h, current_w = current_shape
    
    # Center crop or pad as needed
    if current_h >= target_h and current_w >= target_w:
        # Crop
        start_h = (current_h - target_h) // 2
        start_w = (current_w - target_w) // 2
        return data[start_h:start_h + target_h, start_w:start_w + target_w]
    else:
        # Pad with zeros
        result = np.zeros(target_shape, dtype=data.dtype)
        start_h = (target_h - current_h) // 2
        start_w = (target_w - current_w) // 2
        end_h = start_h + current_h
        end_w = start_w + current_w
        result[start_h:end_h, start_w:end_w] = data
        return result

def create_trilogy_manual_rgb(red_data, green_data, blue_data, output_path,
                             satpercent=0.001, colorsatfac=1,
                             asinh_softening=0.05, brightness_boost=3.0,
                             r_asinh_softening=None, g_asinh_softening=None, b_asinh_softening=None,
                             r_brightness_boost=None, g_brightness_boost=None, b_brightness_boost=None,
                             r_noiselum=None, g_noiselum=None, b_noiselum=None,
                             noiselum=0.15):
    """Create RGB image using TRILOGY-style algorithms with per-channel asinh stretch."""
    if not PILLOW_AVAILABLE:
        print(f"    Error: Pillow library not available. Please install: pip install Pillow")
        return False
    
    try:
        from PIL import Image
        
        # Set per-channel defaults if not specified
        channel_params = {
            'red': {
                'asinh_softening': r_asinh_softening if r_asinh_softening is not None else asinh_softening,
                'brightness_boost': r_brightness_boost if r_brightness_boost is not None else brightness_boost,
                'noiselum': r_noiselum if r_noiselum is not None else noiselum
            },
            'green': {
                'asinh_softening': g_asinh_softening if g_asinh_softening is not None else asinh_softening,
                'brightness_boost': g_brightness_boost if g_brightness_boost is not None else brightness_boost,
                'noiselum': g_noiselum if g_noiselum is not None else noiselum
            },
            'blue': {
                'asinh_softening': b_asinh_softening if b_asinh_softening is not None else asinh_softening,
                'brightness_boost': b_brightness_boost if b_brightness_boost is not None else brightness_boost,
                'noiselum': b_noiselum if b_noiselum is not None else noiselum
            }
        }
        
        # Enhanced scaling algorithm with per-channel asinh stretch
        def trilogy_asinh_scale(data, channel_name, satpercent, params):
            """Apply asinh stretch with per-channel parameters."""
            # Remove invalid data
            valid_data = data[np.isfinite(data)]
            if len(valid_data) == 0:
                return np.zeros_like(data)
            
            # Get channel-specific parameters
            softening = params['asinh_softening']
            brightness = params['brightness_boost']
            channel_noiselum = params['noiselum']
            
            # Calculate robust statistics
            sorted_data = np.sort(valid_data)
            n_pixels = len(sorted_data)
            
            # Calculate noise level (for background subtraction)
            noise_idx = int(channel_noiselum * n_pixels)
            if noise_idx >= n_pixels:
                noise_idx = n_pixels - 1
            noise_level = sorted_data[noise_idx]
            
            # Calculate saturation level for normalization
            sat_idx = int((1.0 - satpercent) * n_pixels)
            if sat_idx >= n_pixels:
                sat_idx = n_pixels - 1
            sat_level = sorted_data[sat_idx]
            
            print(f"      {channel_name.capitalize()} channel: noise={noise_level:.6f}, sat={sat_level:.6f}, "
                  f"softening={softening:.3f}, brightness={brightness:.1f}")
            
            # Background subtraction
            data_bg_sub = data - noise_level
            data_bg_sub = np.maximum(data_bg_sub, 0)
            
            # Normalize to make asinh work properly
            if sat_level > noise_level:
                data_norm = data_bg_sub / (sat_level - noise_level)
            else:
                data_norm = data_bg_sub
            
            # Apply asinh stretch with per-channel brightness boost
            stretched = np.arcsinh(data_norm * brightness / softening) / np.arcsinh(brightness / softening)
            stretched = np.clip(stretched, 0, 1)
            
            return stretched
        
        # Scale each channel using per-channel asinh stretch parameters
        print(f"    Scaling channels with per-channel asinh stretch...")
        red_scaled = trilogy_asinh_scale(red_data, 'red', satpercent, channel_params['red'])
        green_scaled = trilogy_asinh_scale(green_data, 'green', satpercent, channel_params['green'])
        blue_scaled = trilogy_asinh_scale(blue_data, 'blue', satpercent, channel_params['blue'])
        
        # Apply color saturation boost
        if colorsatfac > 1:
            print(f"    Applying color saturation factor: {colorsatfac}")
            luminance = 0.299 * red_scaled + 0.587 * green_scaled + 0.114 * blue_scaled
            
            red_scaled = luminance + colorsatfac * (red_scaled - luminance)
            green_scaled = luminance + colorsatfac * (green_scaled - luminance)
            blue_scaled = luminance + colorsatfac * (blue_scaled - luminance)
            
            # Clip to valid range
            red_scaled = np.clip(red_scaled, 0, 1)
            green_scaled = np.clip(green_scaled, 0, 1)
            blue_scaled = np.clip(blue_scaled, 0, 1)
        
        # Convert to 8-bit RGB
        print(f"    Converting to RGB image...")
        red_8bit = (red_scaled * 255).astype(np.uint8)
        green_8bit = (green_scaled * 255).astype(np.uint8)
        blue_8bit = (blue_scaled * 255).astype(np.uint8)
        
        # Create RGB image (flip vertically to match astronomical convention)
        height, width = red_8bit.shape
        rgb_array = np.zeros((height, width, 3), dtype=np.uint8)
        rgb_array[:, :, 0] = np.flipud(red_8bit)
        rgb_array[:, :, 1] = np.flipud(green_8bit) 
        rgb_array[:, :, 2] = np.flipud(blue_8bit)
        
        # Save as PNG
        image = Image.fromarray(rgb_array)
        image.save(output_path)
        print(f"    RGB image saved: {os.path.basename(output_path)}")
        
        return True
        
    except Exception as e:
        print(f"    Error in RGB creation: {e}")
        return False

def identify_filter_files(fits_files, red_filter='f322w2', green_filter='f322w2', blue_filter='f150w2'):
    """Identify which FITS files correspond to which RGB channels based on filter names (case-insensitive)."""
    red_file = None
    green_file = None  
    blue_file = None
    
    # Convert filter names to lowercase for case-insensitive matching
    red_filter_lower = red_filter.lower()
    green_filter_lower = green_filter.lower()
    blue_filter_lower = blue_filter.lower()
    
    print(f"    Looking for filters: Red='{red_filter}', Green='{green_filter}', Blue='{blue_filter}'")
    
    # Track potential duplicates
    red_matches = []
    green_matches = []
    blue_matches = []
    
    for fits_file in fits_files:
        filename_lower = fits_file.name.lower()
        
        # Check for red filter (case-insensitive)
        if red_filter_lower in filename_lower:
            red_matches.append(fits_file)
        
        # Check for green filter (case-insensitive)  
        if green_filter_lower in filename_lower:
            green_matches.append(fits_file)
            
        # Check for blue filter (case-insensitive)
        if blue_filter_lower in filename_lower:
            blue_matches.append(fits_file)
    
    # Select files and warn about duplicates
    if red_matches:
        red_file = red_matches[0]
        print(f"      Found red channel: {red_file.name}")
        if len(red_matches) > 1:
            print(f"        Warning: Multiple red filter matches found, using first: {[f.name for f in red_matches]}")
    
    if green_matches:
        green_file = green_matches[0]
        print(f"      Found green channel: {green_file.name}")
        if len(green_matches) > 1:
            print(f"        Warning: Multiple green filter matches found, using first: {[f.name for f in green_matches]}")
    
    if blue_matches:
        blue_file = blue_matches[0]
        print(f"      Found blue channel: {blue_file.name}")
        if len(blue_matches) > 1:
            print(f"        Warning: Multiple blue filter matches found, using first: {[f.name for f in blue_matches]}")
    
    return red_file, green_file, blue_file

def create_rgb_from_fits(red_file, green_file, blue_file, output_path, rgb_params):
    """Create RGB image from three FITS files."""
    print(f"  Creating RGB image...")
    print(f"    Red: {red_file.name}")
    print(f"    Green: {green_file.name}")
    print(f"    Blue: {blue_file.name}")
    
    try:
        # Load all files with WCS
        data_list = []
        wcs_list = []
        data_shapes = []
        pixel_scales = []
        
        for fits_file in [red_file, green_file, blue_file]:
            data, header, wcs, scale = load_fits_with_wcs(fits_file)
            data_list.append(data)
            wcs_list.append(wcs)
            data_shapes.append(data.shape)
            pixel_scales.append(scale)
        
        # Find common WCS grid
        ref_wcs, ref_shape = find_common_wcs(wcs_list, data_shapes, pixel_scales)
        
        # Reproject all images to common grid
        print(f"    Reprojecting to common grid...")
        aligned_data = reproject_to_common_grid(data_list, wcs_list, data_shapes, pixel_scales, ref_wcs, ref_shape)
        
        # Apply individual channel scaling if specified
        for i, scale in enumerate([rgb_params.get('r_scale', 1.0), rgb_params.get('g_scale', 1.0), rgb_params.get('b_scale', 1.0)]):
            if scale != 1.0:
                aligned_data[i] = aligned_data[i] * scale
        
        # Create RGB image using manual TRILOGY implementation
        success = create_trilogy_manual_rgb(
            aligned_data[0], aligned_data[1], aligned_data[2], output_path,
            satpercent=rgb_params.get('satpercent', 0.001),
            colorsatfac=rgb_params.get('colorsatfac', 1),
            asinh_softening=rgb_params.get('asinh_softening', 0.05),
            brightness_boost=rgb_params.get('brightness_boost', 3.0),
            r_asinh_softening=rgb_params.get('r_asinh_softening'),
            g_asinh_softening=rgb_params.get('g_asinh_softening'),
            b_asinh_softening=rgb_params.get('b_asinh_softening'),
            r_brightness_boost=rgb_params.get('r_brightness_boost'),
            g_brightness_boost=rgb_params.get('g_brightness_boost'),
            b_brightness_boost=rgb_params.get('b_brightness_boost'),
            r_noiselum=rgb_params.get('r_noiselum'),
            g_noiselum=rgb_params.get('g_noiselum'),
            b_noiselum=rgb_params.get('b_noiselum'),
            noiselum=rgb_params.get('noiselum', 0.15)
        )
        
        return success
        
    except Exception as e:
        print(f"    Error creating RGB: {e}")
        return False

def process_galaxy_directory(input_dir, output_dir, preprocess_params, rgb_params):
    """Process a single galaxy directory (preprocess FITS + create RGB)."""
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    print(f"\nProcessing galaxy directory: {input_path.name}")
    print(f"  Using RGB filter mapping: R='{rgb_params.get('r', 'f322w2')}', G='{rgb_params.get('g', 'f322w2')}', B='{rgb_params.get('b', 'f150w2')}'")
    
    # Find all FITS files
    fits_files = list(input_path.glob('*.fits')) + list(input_path.glob('*.fit'))
    
    if not fits_files:
        print(f"  No FITS files found in {input_dir}")
        return False
    
    print(f"  Found {len(fits_files)} FITS files")
    
    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Preprocess all FITS files
    print(f"  Step 1: Preprocessing FITS files...")
    preprocessed_files = []
    
    for fits_file in fits_files:
        output_fits = output_path / fits_file.name
        success = preprocess_jwst_fits(
            str(fits_file), str(output_fits),
            contrast=preprocess_params.get('contrast', 1.02),
            bias=preprocess_params.get('bias', 0.43),
            max_percent=preprocess_params.get('max_percent', 99)
        )
        
        if success:
            preprocessed_files.append(output_fits)
    
    if not preprocessed_files:
        print(f"  Error: No files were successfully preprocessed")
        return False
    
    # Step 2: Identify RGB channel files
    print(f"  Step 2: Identifying RGB channel files...")
    red_file, green_file, blue_file = identify_filter_files(
        preprocessed_files,
        red_filter=rgb_params.get('r', 'f322w2'),
        green_filter=rgb_params.get('g', 'f322w2'),
        blue_filter=rgb_params.get('b', 'f150w2')
    )
    
    if not all([red_file, green_file, blue_file]):
        print(f"  Warning: Could not identify all RGB channels")
        print(f"    Searched for:")
        print(f"      Red filter:   '{rgb_params.get('r', 'f322w2')}' → {red_file.name if red_file else 'NOT FOUND'}")
        print(f"      Green filter: '{rgb_params.get('g', 'f322w2')}' → {green_file.name if green_file else 'NOT FOUND'}")
        print(f"      Blue filter:  '{rgb_params.get('b', 'f150w2')}' → {blue_file.name if blue_file else 'NOT FOUND'}")
        print(f"    Available files: {[f.name for f in preprocessed_files]}")
        print(f"    Tip: Use -r, -g, -b arguments to specify filter codes (case-insensitive)")
        return False
    
    # Step 3: Create RGB image
    print(f"  Step 3: Creating RGB image...")
    rgb_output = output_path / f"{input_path.name}_RGB.png"
    
    success = create_rgb_from_fits(red_file, green_file, blue_file, str(rgb_output), rgb_params)
    
    if success:
        print(f"  ✓ Successfully processed {input_path.name}")
        return True
    else:
        print(f"  ✗ Failed to create RGB for {input_path.name}")
        return False

def process_all_galaxies(input_dir, output_dir, preprocess_params, rgb_params):
    """Process all galaxy directories in the input directory."""
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    print(f"=== JWST FITS Processing Pipeline ===")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"RGB Filter mapping: Red='{rgb_params.get('r', 'f322w2')}', Green='{rgb_params.get('g', 'f322w2')}', Blue='{rgb_params.get('b', 'f150w2')}' (case-insensitive)")
    print(f"Preprocessing: contrast={preprocess_params.get('contrast', 1.02)}, bias={preprocess_params.get('bias', 0.43)}")
    print(f"RGB stretch: asinh_softening={rgb_params.get('asinh_softening', 0.05)}, brightness_boost={rgb_params.get('brightness_boost', 3.0)}")
    print()
    
    if not input_path.exists():
        print(f"Error: Input directory '{input_dir}' does not exist")
        return
    
    # Find all subdirectories
    galaxy_dirs = [d for d in input_path.iterdir() if d.is_dir()]
    
    if not galaxy_dirs:
        print(f"No subdirectories found in {input_dir}")
        return
    
    print(f"Found {len(galaxy_dirs)} galaxy directories to process:")
    for galaxy_dir in galaxy_dirs:
        print(f"  - {galaxy_dir.name}")
    print()
    
    # Process each galaxy directory
    successful = 0
    failed = 0
    
    for galaxy_dir in galaxy_dirs:
        galaxy_output_dir = output_path / galaxy_dir.name
        
        try:
            success = process_galaxy_directory(galaxy_dir, galaxy_output_dir, preprocess_params, rgb_params)
            if success:
                successful += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  Error processing {galaxy_dir.name}: {e}")
            failed += 1
    
    print(f"\n=== Processing Summary ===")
    print(f"Successfully processed: {successful} galaxy directories")
    print(f"Failed: {failed} galaxy directories")
    print(f"Output directory: {output_dir}")

def main():
    """Main function with command line interface."""
    parser = argparse.ArgumentParser(
        description="Combined JWST FITS preprocessing and RGB image creation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all galaxy directories with default settings (f322w2 for red/green, f150w2 for blue)
  python combined_fits_rgb_processor.py --input-dir ./data --output-dir ./processed_fits_rgb

  # Custom filter mapping using short arguments
  python combined_fits_rgb_processor.py --input-dir ./data --output-dir ./processed_fits_rgb -r f322w2 -g f322w2 -b f150w2

  # Different filter combinations (case-insensitive)
  python combined_fits_rgb_processor.py --input-dir ./data --output-dir ./processed_fits_rgb -r F322W2 -g F150W2 -b F115W

  # Custom RGB parameters for blue channel adjustment
  python combined_fits_rgb_processor.py --input-dir ./data --output-dir ./processed_fits_rgb -r f322w2 -g f322w2 -b f150w2 \\
                                       --b-asinh-softening 0.1 --b-brightness-boost 1.5 --b-noiselum 0.3

  # Conservative preprocessing with custom filters
  python combined_fits_rgb_processor.py --input-dir ./data --output-dir ./processed_fits_rgb -r f444w -g f277w -b f150w \\
                                       --contrast 1.5 --bias 0.5 --max-percent 95

  # High contrast RGB with multiple filter options
  python combined_fits_rgb_processor.py --input-dir ./data --output-dir ./processed_fits_rgb -r F322W2 -g f200w -b F115W \\
                                       --asinh-softening 0.02 --brightness-boost 5.0
        """
    )
    
    # Input/Output directories
    parser.add_argument('--input-dir', required=True, help='Input directory containing galaxy subdirectories')
    parser.add_argument('--output-dir', required=True, help='Output directory for processed files and RGB images')
    
    # RGB Channel filter identification (case-insensitive)
    parser.add_argument('-r', '--red', default='f322w2', help='Filter code for red channel (case-insensitive, default: f322w2)')
    parser.add_argument('-g', '--green', default='f322w2', help='Filter code for green channel (case-insensitive, default: f322w2)')
    parser.add_argument('-b', '--blue', default='f150w2', help='Filter code for blue channel (case-insensitive, default: f150w2)')
    
    # FITS preprocessing parameters
    parser.add_argument('--contrast', type=float, default=1.02, help='DS9 contrast parameter (default: 1.02)')
    parser.add_argument('--bias', type=float, default=0.43, help='DS9 bias parameter (default: 0.43)')
    parser.add_argument('--max-percent', type=float, default=99, help='Maximum percentile for data clipping (default: 99)')
    
    # Channel scaling
    parser.add_argument('--r-scale', type=float, default=1.0, help='Red channel scaling (default: 1.0)')
    parser.add_argument('--g-scale', type=float, default=1.0, help='Green channel scaling (default: 1.0)')
    parser.add_argument('--b-scale', type=float, default=1.0, help='Blue channel scaling (default: 1.0)')
    
    # RGB creation parameters
    parser.add_argument('--satpercent', type=float, default=0.001, help='Percentage of pixels to saturate (default: 0.001)')
    parser.add_argument('--colorsatfac', type=float, default=1, help='Color saturation factor (default: 1)')
    parser.add_argument('--noiselum', type=float, default=0.15, help='Noise luminosity threshold (default: 0.15)')
    
    # Asinh stretch parameters
    parser.add_argument('--asinh-softening', type=float, default=0.05, help='Default asinh softening parameter (default: 0.05)')
    parser.add_argument('--brightness-boost', type=float, default=3.0, help='Default brightness multiplier (default: 3.0)')
    
    # Per-channel asinh stretch parameters
    parser.add_argument('--r-asinh-softening', type=float, default=None, help='Red channel asinh softening')
    parser.add_argument('--g-asinh-softening', type=float, default=None, help='Green channel asinh softening')
    parser.add_argument('--b-asinh-softening', type=float, default=None, help='Blue channel asinh softening')
    
    parser.add_argument('--r-brightness-boost', type=float, default=None, help='Red channel brightness boost')
    parser.add_argument('--g-brightness-boost', type=float, default=None, help='Green channel brightness boost')
    parser.add_argument('--b-brightness-boost', type=float, default=None, help='Blue channel brightness boost')
    
    parser.add_argument('--r-noiselum', type=float, default=None, help='Red channel noise luminosity threshold')
    parser.add_argument('--g-noiselum', type=float, default=None, help='Green channel noise luminosity threshold')
    parser.add_argument('--b-noiselum', type=float, default=None, help='Blue channel noise luminosity threshold')
    
    args = parser.parse_args()
    
    # Suppress warnings
    warnings.filterwarnings('ignore')
    
    # Prepare parameter dictionaries
    preprocess_params = {
        'contrast': args.contrast,
        'bias': args.bias,
        'max_percent': args.max_percent
    }
    
    rgb_params = {
        'r': args.red,
        'g': args.green,
        'b': args.blue,
        'r_scale': args.r_scale,
        'g_scale': args.g_scale,
        'b_scale': args.b_scale,
        'satpercent': args.satpercent,
        'colorsatfac': args.colorsatfac,
        'noiselum': args.noiselum,
        'asinh_softening': args.asinh_softening,
        'brightness_boost': args.brightness_boost,
        'r_asinh_softening': args.r_asinh_softening,
        'g_asinh_softening': args.g_asinh_softening,
        'b_asinh_softening': args.b_asinh_softening,
        'r_brightness_boost': args.r_brightness_boost,
        'g_brightness_boost': args.g_brightness_boost,
        'b_brightness_boost': args.b_brightness_boost,
        'r_noiselum': args.r_noiselum,
        'g_noiselum': args.g_noiselum,
        'b_noiselum': args.b_noiselum
    }
    
    # Process all galaxy directories
    process_all_galaxies(args.input_dir, args.output_dir, preprocess_params, rgb_params)

if __name__ == "__main__":
    main()
