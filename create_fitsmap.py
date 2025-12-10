#!/usr/bin/env python3
"""
FITSMap Generator for Preprocessed JWST Data and RGB Images

This script creates interactive FITSMap tiles from preprocessed FITS files
and RGB images (PNG/JPG), preserving WCS coordinates for proper zooming 
and catalog integration.

Usage:
    python create_fitsmap.py input_directory output_directory
    python create_fitsmap.py input_directory output_directory --serve
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
from fitsmap import convert
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import numpy as np
import tempfile
import re

def get_wcs_and_shape(fits_file):
    """Get WCS and data shape from a FITS file."""
    try:
        with fits.open(fits_file) as hdul:
            # Find extension with data
            for i, hdu in enumerate(hdul):
                if hdu.data is not None:
                    header = hdu.header.copy()
                    shape = hdu.data.shape
                    try:
                        wcs = WCS(header)
                        return wcs, shape, i
                    except:
                        continue
        return None, None, None
    except:
        return None, None, None

def find_reference_wcs(fits_files):
    """Find the FITS file with finest pixel scale to use as reference."""
    print("Analyzing WCS information for alignment...")
    
    best_file = None
    finest_scale = float('inf')
    best_wcs = None
    best_shape = None
    
    for fits_file in fits_files:
        wcs, shape, ext = get_wcs_and_shape(fits_file)
        if wcs is not None:
            try:
                scales = wcs.proj_plane_pixel_scales()
                if hasattr(scales[0], 'unit'):
                    if scales[0].unit.to_string() == 'deg':
                        scale_arcsec = scales[0].to_value('arcsec')
                    else:
                        scale_arcsec = scales[0].value * 3600
                else:
                    scale_arcsec = float(scales[0]) * 3600
                
                # Also check WCS center coordinates for debugging
                center_ra = wcs.wcs.crval[0] if hasattr(wcs.wcs, 'crval') else 'unknown'
                center_dec = wcs.wcs.crval[1] if hasattr(wcs.wcs, 'crval') else 'unknown'
                
                print(f"  {os.path.basename(fits_file)}: {scale_arcsec:.3f} arcsec/pixel, shape {shape}")
                print(f"    Center: RA={center_ra:.6f}°, Dec={center_dec:.6f}°")
                
                if scale_arcsec < finest_scale:
                    finest_scale = scale_arcsec
                    best_file = fits_file
                    best_wcs = wcs
                    best_shape = shape
                    
            except Exception as e:
                print(f"  Warning: Could not get pixel scale for {fits_file}: {e}")
    
    if best_file:
        print(f"Using {os.path.basename(best_file)} as reference ({finest_scale:.3f} arcsec/pixel)")
    else:
        print("Warning: No suitable reference file found!")
    
    return best_file, best_wcs, best_shape

def add_coordinate_search(output_dir, reference_file):
    """
    Add coordinate search functionality to the generated FITSMap.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing the generated FITSMap
    reference_file : str
        Reference FITS file for WCS information
    """
    index_html = os.path.join(output_dir, 'index.html')
    
    if not os.path.exists(index_html):
        print("Warning: index.html not found, skipping coordinate search addition")
        return
    
    print("Adding coordinate search functionality...")
    
    # Extract WCS parameters from reference file
    try:
        with fits.open(reference_file) as hdul:
            # Find extension with data
            wcs_header = None
            for hdu in hdul:
                if hdu.data is not None:
                    wcs_header = hdu.header
                    break
            
            if wcs_header is None:
                print("Warning: Could not extract WCS from reference file")
                return
            
            wcs = WCS(wcs_header)
            
            # Get key WCS parameters for JavaScript
            crval1 = wcs_header.get('CRVAL1', 0)
            crval2 = wcs_header.get('CRVAL2', 0)
            crpix1 = wcs_header.get('CRPIX1', 0)
            crpix2 = wcs_header.get('CRPIX2', 0)
            
            # Get CD matrix or CDELT values
            cd1_1 = wcs_header.get('CD1_1', wcs_header.get('CDELT1', 0))
            cd1_2 = wcs_header.get('CD1_2', 0)
            cd2_1 = wcs_header.get('CD2_1', 0)
            cd2_2 = wcs_header.get('CD2_2', wcs_header.get('CDELT2', 0))
            
    except Exception as e:
        print(f"Warning: Could not extract WCS parameters from {reference_file}: {e}")
        return
    
    # Read the existing HTML
    try:
        with open(index_html, 'r', encoding='utf-8') as f:
            html_content = f.read()
    except Exception as e:
        print(f"Warning: Could not read index.html: {e}")
        return
    
    # HTML for coordinate search box
    coordinate_search_html = f"""
    <!-- Coordinate Search Box -->
    <div id="coordinate-search" style="
        position: absolute; 
        bottom: 10px; 
        right: 10px; 
        background: rgba(255, 255, 255, 0.92); 
        padding: 6px; 
        border-radius: 4px; 
        box-shadow: 0 2px 4px rgba(0,0,0,0.3);
        font-family: Arial, sans-serif;
        font-size: 11px;
        z-index: 1000;
        width: 160px;">
        <div style="margin-bottom: 4px; font-weight: bold; font-size: 10px;">Go to Coordinates:</div>
        <div style="margin-bottom: 4px;">
            <input type="text" id="ra-input" placeholder="RA (deg or hh:mm:ss)" 
                   style="width: 100%; padding: 2px; margin-bottom: 2px; border: 1px solid #ccc; border-radius: 2px; font-size: 10px; box-sizing: border-box;">
            <input type="text" id="dec-input" placeholder="Dec (deg or dd:mm:ss)" 
                   style="width: 100%; padding: 2px; margin-bottom: 2px; border: 1px solid #ccc; border-radius: 2px; font-size: 10px; box-sizing: border-box;">
        </div>
        <button onclick="goToCoordinates()" style="
            width: 100%; 
            padding: 3px; 
            background: #007cba; 
            color: white; 
            border: none; 
            border-radius: 2px; 
            cursor: pointer;
            font-size: 10px;">
            Go
        </button>
        <div style="margin-top: 3px; font-size: 8px; color: #666; line-height: 1.2;">
            Examples: 150.5, -20.3<br>
            or 10:02:00, -20:18:00
        </div>
    </div>
    """
    
def add_coordinate_search(output_dir, reference_file):
    """
    Add coordinate search functionality to the generated FITSMap.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing the generated FITSMap
    reference_file : str
        Reference FITS file for WCS information
    """
    index_html = os.path.join(output_dir, 'index.html')
    
    if not os.path.exists(index_html):
        print("Warning: index.html not found, skipping coordinate search addition")
        return
    
    print("Adding coordinate search functionality...")
    
    # Extract WCS parameters from reference file and generate sample coordinate grid
    try:
        with fits.open(reference_file) as hdul:
            # Find extension with data
            wcs_header = None
            for hdu in hdul:
                if hdu.data is not None:
                    wcs_header = hdu.header
                    data_shape = hdu.data.shape
                    break
            
            if wcs_header is None:
                print("Warning: Could not extract WCS from reference file")
                return
            
            wcs = WCS(wcs_header)
            
            # Generate a grid of pixel coordinates to create a lookup table
            # This provides more accurate transformation than manual calculation
            print("Generating coordinate transformation grid...")
            
            # Create a grid spanning the image
            nx, ny = data_shape[1], data_shape[0]  # astropy uses (ny, nx) order
            
            # Create sample points across the image for transformation lookup
            grid_points = 20  # 20x20 grid for interpolation
            x_samples = np.linspace(0, nx-1, grid_points)
            y_samples = np.linspace(0, ny-1, grid_points)
            xx, yy = np.meshgrid(x_samples, y_samples)
            
            # Convert to world coordinates using proper astropy transformation
            # Use 0-based indexing for astropy
            world_coords = wcs.pixel_to_world(xx.flatten(), yy.flatten())
            
            # Extract RA/Dec in degrees
            ra_grid = world_coords.ra.degree.reshape(grid_points, grid_points)
            dec_grid = world_coords.dec.degree.reshape(grid_points, grid_points)
            
            # Get basic WCS parameters for fallback
            crval1 = wcs_header.get('CRVAL1', 0)
            crval2 = wcs_header.get('CRVAL2', 0)
            crpix1 = wcs_header.get('CRPIX1', 0)
            crpix2 = wcs_header.get('CRPIX2', 0)
            
            # Get CD matrix or CDELT values
            cd1_1 = wcs_header.get('CD1_1', wcs_header.get('CDELT1', 0))
            cd1_2 = wcs_header.get('CD1_2', 0)
            cd2_1 = wcs_header.get('CD2_1', 0)
            cd2_2 = wcs_header.get('CD2_2', wcs_header.get('CDELT2', 0))
            
            print(f"Image shape: {nx} x {ny} pixels")
            print(f"WCS center: RA={crval1:.6f}°, Dec={crval2:.6f}°")
            print(f"Generated {grid_points}x{grid_points} coordinate grid for accurate transformation")
            
    except Exception as e:
        print(f"Warning: Could not extract WCS parameters from {reference_file}: {e}")
        return
    
    # Read the existing HTML
    try:
        with open(index_html, 'r', encoding='utf-8') as f:
            html_content = f.read()
    except Exception as e:
        print(f"Warning: Could not read index.html: {e}")
        return
    
    # HTML for coordinate search box
    coordinate_search_html = f"""
    <!-- Coordinate Search Box -->
    <div id="coordinate-search" style="
        position: absolute; 
        bottom: 10px; 
        right: 10px; 
        background: rgba(255, 255, 255, 0.92); 
        padding: 6px; 
        border-radius: 4px; 
        box-shadow: 0 2px 4px rgba(0,0,0,0.3);
        font-family: Arial, sans-serif;
        font-size: 11px;
        z-index: 1000;
        width: 160px;">
        <div style="margin-bottom: 4px; font-weight: bold; font-size: 10px;">Go to Coordinates:</div>
        <div style="margin-bottom: 4px;">
            <input type="text" id="ra-input" placeholder="RA (deg or hh:mm:ss)" 
                   style="width: 100%; padding: 2px; margin-bottom: 2px; border: 1px solid #ccc; border-radius: 2px; font-size: 10px; box-sizing: border-box;">
            <input type="text" id="dec-input" placeholder="Dec (deg or dd:mm:ss)" 
                   style="width: 100%; padding: 2px; margin-bottom: 2px; border: 1px solid #ccc; border-radius: 2px; font-size: 10px; box-sizing: border-box;">
        </div>
        <button onclick="goToCoordinates()" style="
            width: 100%; 
            padding: 3px; 
            background: #007cba; 
            color: white; 
            border: none; 
            border-radius: 2px; 
            cursor: pointer;
            font-size: 10px;">
            Go
        </button>
        <div style="margin-top: 3px; font-size: 8px; color: #666; line-height: 1.2;">
            Examples: 164.38, 58.00<br>
            or 10:57:30, +58:00:00
        </div>
    </div>
    """
    
    # Convert numpy arrays to JavaScript arrays
    x_samples_js = [float(x) for x in x_samples]
    y_samples_js = [float(y) for y in y_samples]
    ra_grid_js = [[float(val) for val in row] for row in ra_grid]
    dec_grid_js = [[float(val) for val in row] for row in dec_grid]
    
def add_coordinate_search(output_dir, reference_file):
    """
    Add coordinate search functionality to the generated FITSMap.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing the generated FITSMap
    reference_file : str
        Reference FITS file for WCS information
    """
    index_html = os.path.join(output_dir, 'index.html')
    
    if not os.path.exists(index_html):
        print("Warning: index.html not found, skipping coordinate search addition")
        return
    
    print("Adding coordinate search functionality...")
    
    # Extract basic WCS parameters from reference file
    try:
        with fits.open(reference_file) as hdul:
            # Find extension with data
            wcs_header = None
            for hdu in hdul:
                if hdu.data is not None:
                    wcs_header = hdu.header
                    data_shape = hdu.data.shape
                    break
            
            if wcs_header is None:
                print("Warning: Could not extract WCS from reference file")
                return
            
            wcs = WCS(wcs_header)
            
            # Get image center coordinates as a reference point
            center_x = data_shape[1] / 2.0
            center_y = data_shape[0] / 2.0
            center_coord = wcs.pixel_to_world(center_x, center_y)
            center_ra = center_coord.ra.degree
            center_dec = center_coord.dec.degree
            
            # Get pixel scale (degrees per pixel)
            pixel_scales = wcs.proj_plane_pixel_scales()
            if hasattr(pixel_scales[0], 'value'):
                pixel_scale_ra = abs(pixel_scales[0].value)
                pixel_scale_dec = abs(pixel_scales[1].value)
            else:
                pixel_scale_ra = abs(pixel_scales[0])
                pixel_scale_dec = abs(pixel_scales[1])
            
            # Get basic WCS parameters
            crval1 = wcs_header.get('CRVAL1', center_ra)
            crval2 = wcs_header.get('CRVAL2', center_dec)
            crpix1 = wcs_header.get('CRPIX1', center_x)
            crpix2 = wcs_header.get('CRPIX2', center_y)
            
            print(f"Image shape: {data_shape[1]} x {data_shape[0]} pixels")
            print(f"Image center: RA={center_ra:.6f}°, Dec={center_dec:.6f}°")
            print(f"Pixel scale: {pixel_scale_ra*3600:.3f}\" x {pixel_scale_dec*3600:.3f}\" per pixel")
            
    except Exception as e:
        print(f"Warning: Could not extract WCS parameters from {reference_file}: {e}")
        return
    
    # Read the existing HTML
    try:
        with open(index_html, 'r', encoding='utf-8') as f:
            html_content = f.read()
    except Exception as e:
        print(f"Warning: Could not read index.html: {e}")
        return
    
    # HTML for coordinate search box
    coordinate_search_html = f"""
    <!-- Coordinate Search Box -->
    <div id="coordinate-search" style="
        position: absolute; 
        bottom: 10px; 
        right: 10px; 
        background: rgba(255, 255, 255, 0.92); 
        padding: 6px; 
        border-radius: 4px; 
        box-shadow: 0 2px 4px rgba(0,0,0,0.3);
        font-family: Arial, sans-serif;
        font-size: 11px;
        z-index: 1000;
        width: 160px;">
        <div style="margin-bottom: 4px; font-weight: bold; font-size: 10px;">Go to Coordinates:</div>
        <div style="margin-bottom: 4px;">
            <input type="text" id="ra-input" placeholder="RA (deg or hh:mm:ss)" 
                   style="width: 100%; padding: 2px; margin-bottom: 2px; border: 1px solid #ccc; border-radius: 2px; font-size: 10px; box-sizing: border-box;"
                   value="{center_ra:.6f}">
            <input type="text" id="dec-input" placeholder="Dec (deg or dd:mm:ss)" 
                   style="width: 100%; padding: 2px; margin-bottom: 2px; border: 1px solid #ccc; border-radius: 2px; font-size: 10px; box-sizing: border-box;"
                   value="{center_dec:.6f}">
        </div>
        <button onclick="goToCoordinates()" style="
            width: 100%; 
            padding: 3px; 
            background: #007cba; 
            color: white; 
            border: none; 
            border-radius: 2px; 
            cursor: pointer;
            font-size: 10px;">
            Go
        </button>
        <div style="margin-top: 3px; font-size: 8px; color: #666; line-height: 1.2;">
            Click to test with center coords<br>
            or enter your own RA/Dec
        </div>
    </div>
    """
    
def add_coordinate_search(output_dir, reference_file):
    """
    Add coordinate search functionality to the generated FITSMap.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing the generated FITSMap
    reference_file : str
        Reference FITS file for WCS information
    """
    index_html = os.path.join(output_dir, 'index.html')
    
    if not os.path.exists(index_html):
        print("Warning: index.html not found, skipping coordinate search addition")
        return
    
    print("Adding coordinate search functionality...")
    
    # Extract WCS and generate coordinate conversion data
    try:
        with fits.open(reference_file) as hdul:
            # Find extension with data
            wcs_header = None
            for hdu in hdul:
                if hdu.data is not None:
                    wcs_header = hdu.header
                    data_shape = hdu.data.shape
                    break
            
            if wcs_header is None:
                print("Warning: Could not extract WCS from reference file")
                return
            
            wcs = WCS(wcs_header)
            
            # Get image dimensions and center
            nx, ny = data_shape[1], data_shape[0]  # (width, height)
            center_x = nx / 2.0
            center_y = ny / 2.0
            
            # Get center coordinates
            center_coord = wcs.pixel_to_world(center_x, center_y)
            center_ra = center_coord.ra.degree
            center_dec = center_coord.dec.degree
            
            print(f"Image shape: {nx} x {ny} pixels")
            print(f"Image center: RA={center_ra:.6f}°, Dec={center_dec:.6f}°")
            print(f"Pixel center: ({center_x:.1f}, {center_y:.1f})")
            
            # Test coordinate conversion with a few known points
            test_coords = [
                (center_ra, center_dec),
                (center_ra + 0.01, center_dec),  # Offset by ~36 arcsec
                (center_ra - 0.01, center_dec),
                (center_ra, center_dec + 0.01),
                (center_ra, center_dec - 0.01)
            ]
            
            conversion_examples = []
            for ra, dec in test_coords:
                try:
                    # Convert to pixel coordinates (use world_to_pixel for proper conversion)
                    px, py = wcs.world_to_pixel(ra, dec)
                    # Store both regular and Y-flipped coordinates for testing
                    conversion_examples.append([ra, dec, float(px), float(py), float(nx), float(ny)])
                    print(f"  RA={ra:.6f}°, Dec={dec:.6f}° -> pixel ({px:.2f}, {py:.2f})")
                except Exception as e:
                    print(f"  Error converting RA={ra:.6f}°, Dec={dec:.6f}°: {e}")
            
            # Ensure we have at least one example (the center)
            if not conversion_examples:
                print("Warning: All coordinate conversions failed, using center pixel as fallback")
                conversion_examples = [[center_ra, center_dec, center_x, center_y, nx, ny]]
            
    except Exception as e:
        print(f"Warning: Could not extract WCS parameters from {reference_file}: {e}")
        return
    
    # Read the existing HTML
    try:
        with open(index_html, 'r', encoding='utf-8') as f:
            html_content = f.read()
    except Exception as e:
        print(f"Warning: Could not read index.html: {e}")
        return
    
    # HTML for coordinate search box
    coordinate_search_html = f"""
    <!-- Coordinate Search Box -->
    <div id="coordinate-search" style="
        position: absolute; 
        bottom: 10px; 
        right: 10px; 
        background: rgba(255, 255, 255, 0.92); 
        padding: 6px; 
        border-radius: 4px; 
        box-shadow: 0 2px 4px rgba(0,0,0,0.3);
        font-family: Arial, sans-serif;
        font-size: 11px;
        z-index: 1000;
        width: 160px;">
        <div style="margin-bottom: 4px; font-weight: bold; font-size: 10px;">Go to Coordinates:</div>
        <div style="margin-bottom: 4px;">
            <input type="text" id="ra-input" placeholder="RA (deg or hh:mm:ss)" 
                   style="width: 100%; padding: 2px; margin-bottom: 2px; border: 1px solid #ccc; border-radius: 2px; font-size: 10px; box-sizing: border-box;"
                   value="{center_ra:.6f}">
            <input type="text" id="dec-input" placeholder="Dec (deg or dd:mm:ss)" 
                   style="width: 100%; padding: 2px; margin-bottom: 2px; border: 1px solid #ccc; border-radius: 2px; font-size: 10px; box-sizing: border-box;"
                   value="{center_dec:.6f}">
        </div>
        <button onclick="goToCoordinates()" style="
            width: 100%; 
            padding: 3px; 
            background: #007cba; 
            color: white; 
            border: none; 
            border-radius: 2px; 
            cursor: pointer;
            font-size: 10px;">
            Go
        </button>
        <div style="margin-top: 3px; font-size: 8px; color: #666; line-height: 1.2;">
            Click to test center<br>
            Try: RA ± 0.01°, Dec ± 0.01°
        </div>
    </div>
    """
    
    # JavaScript with proper coordinate conversion
    coordinate_search_js = f"""
    <script>
    // Coordinate conversion examples from Python WCS
    const COORDINATE_EXAMPLES = {conversion_examples};
    
    function parseDegreeOrHMS(input, isRA) {{
        input = input.trim();
        
        // If it contains colons, treat as HMS/DMS
        if (input.includes(':')) {{
            const parts = input.split(':').map(p => parseFloat(p));
            if (parts.length >= 2) {{
                let degrees = Math.abs(parts[0]) + parts[1]/60;
                if (parts.length >= 3) {{
                    degrees += parts[2]/3600;
                }}
                // Handle negative declination
                if (input.startsWith('-') || parts[0] < 0) {{
                    degrees = -degrees;
                }}
                // Convert RA hours to degrees
                if (isRA) {{
                    degrees *= 15;
                }}
                return degrees;
            }}
        }}
        
        // Otherwise treat as decimal degrees
        return parseFloat(input);
    }}
    
    function worldToPixel(ra, dec) {{
        // Check if we have valid coordinate examples
        if (!COORDINATE_EXAMPLES || COORDINATE_EXAMPLES.length === 0) {{
            console.warn('No coordinate examples available, using simple fallback');
            return [{center_x}, {center_y}];
        }}
        
        // Use the coordinate examples to interpolate or find nearest match
        let bestDistance = Infinity;
        let bestPixelX = {center_x};
        let bestPixelY = {center_y};
        
        // First, try to find exact matches or very close ones
        for (let example of COORDINATE_EXAMPLES) {{
            if (!example || example.length < 4) {{
                console.warn('Invalid coordinate example:', example);
                continue;
            }}
            
            const [exampleRA, exampleDec, pixelX, pixelY] = example;
            
            const distance = Math.sqrt(
                Math.pow(ra - exampleRA, 2) + Math.pow(dec - exampleDec, 2)
            );
            
            if (distance < bestDistance) {{
                bestDistance = distance;
                bestPixelX = pixelX;
                bestPixelY = pixelY;
            }}
        }}
        
        // If we have a close match (within 0.001 degrees ~ 3.6 arcsec), use it
        if (bestDistance < 0.001) {{
            console.log(`Found close match: distance=${{(bestDistance*3600).toFixed(2)}}" -> pixel (${{bestPixelX.toFixed(2)}}, ${{bestPixelY.toFixed(2)}})`);
            return [bestPixelX, bestPixelY];
        }}
        
        // Otherwise, do linear interpolation from center
        const centerExample = COORDINATE_EXAMPLES[0];
        if (!centerExample || centerExample.length < 4) {{
            console.warn('Invalid center example, using image center');
            return [{center_x}, {center_y}];
        }}
        
        const [centerRA, centerDec, centerPixelX, centerPixelY] = centerExample;
        
        // Find the best scaling using other examples
        let scaleX = -1000.0;  // Default scale (degrees to pixels)
        let scaleY = 1000.0;   // Positive for Dec (usually Y increases with Dec)
        
        for (let i = 1; i < COORDINATE_EXAMPLES.length; i++) {{
            const example = COORDINATE_EXAMPLES[i];
            if (!example || example.length < 4) continue;
            
            const [exampleRA, exampleDec, pixelX, pixelY] = example;
            
            const deltaRA = exampleRA - centerRA;
            const deltaDec = exampleDec - centerDec;
            const deltaPixelX = pixelX - centerPixelX;
            const deltaPixelY = pixelY - centerPixelY;
            
            if (Math.abs(deltaRA) > 1e-6) {{
                scaleX = deltaPixelX / deltaRA;
                console.log(`Scale X from example ${{i}}: ΔRA=${{deltaRA.toFixed(6)}}° -> Δpx=${{deltaPixelX.toFixed(2)}} => scale=${{scaleX.toFixed(1)}}`);
            }}
            if (Math.abs(deltaDec) > 1e-6) {{
                scaleY = deltaPixelY / deltaDec;
                console.log(`Scale Y from example ${{i}}: ΔDec=${{deltaDec.toFixed(6)}}° -> Δpx=${{deltaPixelY.toFixed(2)}} => scale=${{scaleY.toFixed(1)}}`);
            }}
        }}
        
        // Apply linear scaling from center
        const deltaRA = ra - centerRA;
        const deltaDec = dec - centerDec;
        
        const pixelX = centerPixelX + deltaRA * scaleX;
        const pixelY = centerPixelY + deltaDec * scaleY;
        
        console.log(`Linear interpolation: ΔRA=${{deltaRA.toFixed(6)}}°, ΔDec=${{deltaDec.toFixed(6)}}° -> Δpixel(${{(deltaRA*scaleX).toFixed(2)}}, ${{(deltaDec*scaleY).toFixed(2)}})`);
        
        return [pixelX, pixelY];
    }}
    
    function goToCoordinates() {{
        const raInput = document.getElementById('ra-input').value;
        const decInput = document.getElementById('dec-input').value;
        
        if (!raInput || !decInput) {{
            alert('Please enter both RA and Dec coordinates');
            return;
        }}
        
        try {{
            const ra = parseDegreeOrHMS(raInput, true);
            const dec = parseDegreeOrHMS(decInput, false);
            
            if (isNaN(ra) || isNaN(dec)) {{
                alert('Invalid coordinates. Please check your input.');
                return;
            }}
            
            console.log(`\\n=== COORDINATE SEARCH ==`);
            console.log(`Input: RA=${{ra.toFixed(6)}}°, Dec=${{dec.toFixed(6)}}°`);
            
            // Convert world coordinates to pixel coordinates
            const [pixelX, pixelY] = worldToPixel(ra, dec);
            console.log(`Calculated pixel coordinates: (${{pixelX.toFixed(2)}}, ${{pixelY.toFixed(2)}})`);
            
            // Find the map object
            let theMap = null;
            if (typeof map !== 'undefined') {{
                theMap = map;
            }} else {{
                const candidates = ['map', 'leafletMap', 'myMap', 'viewer', 'fitsMap'];
                for (let candidate of candidates) {{
                    if (typeof window[candidate] !== 'undefined' && 
                        window[candidate] && 
                        typeof window[candidate].setView === 'function') {{
                        theMap = window[candidate];
                        window.map = theMap;
                        break;
                    }}
                }}
            }}
            
            if (!theMap) {{
                alert('Cannot find map object. Please wait for page to load completely.');
                return;
            }}
            
            // CRITICAL: FITSMap uses pixel coordinates directly as lat/lng
            // So we use pixelY as lat and pixelX as lng
            const leafletCoords = L.latLng(pixelY, pixelX);
            
            console.log(`Setting Leaflet view to: lat=${{pixelY.toFixed(2)}}, lng=${{pixelX.toFixed(2)}}`);
            
            // Center the map on these pixel coordinates
            theMap.setView(leafletCoords, theMap.getZoom());
            
            // Add a marker to show where we centered
            if (typeof L !== 'undefined') {{
                if (window.searchMarker) {{
                    theMap.removeLayer(window.searchMarker);
                }}
                
                window.searchMarker = L.marker(leafletCoords).addTo(theMap)
                    .bindPopup(`RA: ${{ra.toFixed(6)}}°<br>Dec: ${{dec.toFixed(6)}}°<br>Pixel: (${{pixelX.toFixed(1)}}, ${{pixelY.toFixed(1)}})`)
                    .openPopup();
                
                // Remove marker after 8 seconds
                setTimeout(() => {{
                    if (window.searchMarker) {{
                        theMap.removeLayer(window.searchMarker);
                        window.searchMarker = null;
                    }}
                }}, 8000);
            }}
            
            console.log(`Successfully centered on RA: ${{ra.toFixed(6)}}°, Dec: ${{dec.toFixed(6)}}°`);
            console.log('=== END SEARCH ===\\n');
            
        }} catch (error) {{
            alert('Error: ' + error.message);
            console.error('Error in coordinate search:', error);
        }}
    }}
    
    // Initialize
    document.addEventListener('DOMContentLoaded', function() {{
        const inputs = ['ra-input', 'dec-input'];
        inputs.forEach(id => {{
            const input = document.getElementById(id);
            if (input) {{
                input.addEventListener('keypress', function(e) {{
                    if (e.key === 'Enter') {{
                        goToCoordinates();
                    }}
                }});
            }}
        }});
        
        console.log('FITSMap coordinate search initialized');
        console.log('Image center: RA={center_ra:.6f}°, Dec={center_dec:.6f}°');
        console.log('Coordinate examples available:', COORDINATE_EXAMPLES.length);
        
        // Debug: Show all coordinate examples
        if (COORDINATE_EXAMPLES.length > 0) {{
            console.log('Coordinate conversion examples:');
            COORDINATE_EXAMPLES.forEach((example, i) => {{
                if (example && example.length >= 4) {{
                    const [ra, dec, px, py] = example;
                    console.log(`  ${{i}}: RA=${{ra.toFixed(6)}}°, Dec=${{dec.toFixed(6)}}° -> pixel (${{px.toFixed(2)}}, ${{py.toFixed(2)}})`);
                }} else {{
                    console.log(`  ${{i}}: Invalid example:`, example);
                }}
            }});
        }} else {{
            console.warn('No coordinate examples available - coordinate search may not work properly');
        }}
        
        console.log('Try clicking "Go" to test center coordinates');
    }});
    </script>
    """
    
    # Insert the HTML before the closing body tag
    if '</body>' in html_content:
        html_content = html_content.replace('</body>', coordinate_search_html + '\n</body>')
    else:
        # Fallback: append to end of HTML
        html_content += coordinate_search_html
    
    # Insert the JavaScript before the closing head tag
    if '</head>' in html_content:
        html_content = html_content.replace('</head>', coordinate_search_js + '\n</head>')
    else:
        # Fallback: append to end of HTML
        html_content += coordinate_search_js
    
    # Write the modified HTML back
    try:
        with open(index_html, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print("Coordinate search functionality added successfully!")
        print("You can now enter RA/Dec coordinates in the bottom-left box to center the map")
    except Exception as e:
        print(f"Warning: Could not write modified HTML: {e}")

def find_fits_and_image_files(input_dir):
    """
    Find all FITS and image files in the input directory.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing preprocessed FITS files and RGB images
        
    Returns:
    --------
    tuple
        (fits_files, image_files) - Lists of paths to FITS and image files
    """
    input_path = Path(input_dir)
    
    # Find FITS files
    fits_files = list(input_path.glob('*.fits')) + list(input_path.glob('*.fit'))
    fits_files = [str(f) for f in fits_files]
    
    # Find image files (PNG, JPG, JPEG)
    image_files = (list(input_path.glob('*.png')) + 
                   list(input_path.glob('*.jpg')) + 
                   list(input_path.glob('*.jpeg')))
    image_files = [str(f) for f in image_files]
    
    return fits_files, image_files

def check_wcs_alignment(fits_files):
    """
    Check WCS alignment between FITS files and report any offsets.
    
    Parameters:
    -----------
    fits_files : list
        List of FITS file paths
    """
    print("\nChecking WCS alignment between files:")
    
    wcs_info = []
    for fits_file in fits_files:
        wcs, shape, ext = get_wcs_and_shape(fits_file)
        if wcs is not None and hasattr(wcs.wcs, 'crval'):
            center_ra = wcs.wcs.crval[0]
            center_dec = wcs.wcs.crval[1]
            wcs_info.append({
                'file': fits_file,
                'ra': center_ra,
                'dec': center_dec,
                'wcs': wcs,
                'shape': shape
            })
            print(f"  {os.path.basename(fits_file)}: RA={center_ra:.6f}°, Dec={center_dec:.6f}°")
        else:
            print(f"  {os.path.basename(fits_file)}: No valid WCS found")
    
    if len(wcs_info) < 2:
        print("  Need at least 2 files with valid WCS for alignment check")
        return
    
    # Calculate offsets from first file
    ref_info = wcs_info[0]
    print(f"\nOffsets relative to {os.path.basename(ref_info['file'])}:")
    
    max_offset = 0
    for info in wcs_info[1:]:
        ra_offset = info['ra'] - ref_info['ra']
        dec_offset = info['dec'] - ref_info['dec']
        total_offset_arcsec = np.sqrt(ra_offset**2 + dec_offset**2) * 3600
        max_offset = max(max_offset, total_offset_arcsec)
        
        print(f"  {os.path.basename(info['file'])}: "
              f"ΔRA={ra_offset:.6f}° ({ra_offset*3600:.2f}\"), "
              f"ΔDec={dec_offset:.6f}° ({dec_offset*3600:.2f}\"), "
              f"Total={total_offset_arcsec:.2f}\"")
    
    if max_offset > 1.0:
        print(f"\nSignificant WCS offsets detected (max: {max_offset:.2f}\")")
        print("Recommendation: Enable WCS alignment (default) or manually check your files")
    else:
        print(f"\nSmall WCS offsets (max: {max_offset:.2f}\")")
        print("Files appear to be well-aligned already")
    """
    Select the best FITS file to use as WCS reference.
    Prioritizes files with the most complete WCS information.
    
    Parameters:
    -----------
    fits_files : list
        List of FITS file paths
        
    Returns:
    --------
    str
        Path to the best WCS reference file
    """
    best_file = None
    max_wcs_keywords = 0
    
    wcs_keywords = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CTYPE1', 'CTYPE2']
    
    for fits_file in fits_files:
        try:
            with fits.open(fits_file) as hdul:
                header = hdul[0].header
                wcs_count = sum(1 for keyword in wcs_keywords if keyword in header)
                
                if wcs_count > max_wcs_keywords:
                    max_wcs_keywords = wcs_count
                    best_file = fits_file
                    
        except Exception as e:
            print(f"Warning: Could not read {fits_file}: {e}")
            continue
    
    if best_file is None and fits_files:
        # Fallback to first file if no WCS info found
        best_file = fits_files[0]
        print(f"Warning: No WCS keywords found, using {best_file} as reference")
    
    return best_file

def reproject_fits_file(input_file, reference_wcs, reference_shape, output_file):
    """Reproject a FITS file to match reference WCS and save aligned version."""
    print(f"    Processing {os.path.basename(input_file)}...")
    
    try:
        with fits.open(input_file) as hdul:
            # Find data extension
            data_hdu = None
            for i, hdu in enumerate(hdul):
                if hdu.data is not None:
                    data_hdu = hdu
                    ext_num = i
                    break
            
            if data_hdu is None:
                raise ValueError("No data found in FITS file")
            
            data = data_hdu.data
            header = data_hdu.header.copy()
            
            # Get input WCS and check for issues
            try:
                input_wcs = WCS(header)
                
                # Debug: Print WCS centers to check for offsets
                input_center = input_wcs.wcs.crval if hasattr(input_wcs.wcs, 'crval') else None
                ref_center = reference_wcs.wcs.crval if hasattr(reference_wcs.wcs, 'crval') else None
                
                if input_center is not None and ref_center is not None:
                    ra_offset = input_center[0] - ref_center[0]
                    dec_offset = input_center[1] - ref_center[1]
                    print(f"      Input center: RA={input_center[0]:.6f}°, Dec={input_center[1]:.6f}°")
                    print(f"      Reference center: RA={ref_center[0]:.6f}°, Dec={ref_center[1]:.6f}°")
                    print(f"      Offset: ΔRA={ra_offset:.6f}°, ΔDec={dec_offset:.6f}°")
                    
                    # Check if offset is significant (more than 1 arcsec)
                    total_offset_arcsec = np.sqrt(ra_offset**2 + dec_offset**2) * 3600
                    if total_offset_arcsec > 1.0:
                        print(f"      Significant offset detected: {total_offset_arcsec:.2f} arcsec - reprojecting...")
                    else:
                        print(f"      Small offset: {total_offset_arcsec:.2f} arcsec - reprojecting anyway...")
                
                print(f"      Reprojecting from shape {data.shape} to {reference_shape}...")
                
                # Reproject to reference grid
                reprojected_data, footprint = reproject_interp(
                    (data, input_wcs), reference_wcs, 
                    shape_out=reference_shape, 
                    order='bilinear'
                )
                
                # Check reprojection quality
                valid_pixels = np.sum(np.isfinite(reprojected_data))
                total_pixels = np.prod(reprojected_data.shape)
                valid_fraction = valid_pixels / total_pixels
                
                print(f"      Reprojection completed: {valid_pixels}/{total_pixels} valid pixels ({valid_fraction:.1%})")
                
                if valid_fraction < 0.1:
                    print(f"      Warning: Very few valid pixels after reprojection - may indicate WCS issues")
                
                # Handle NaN values - replace with zeros
                nan_mask = ~np.isfinite(reprojected_data)
                if np.any(nan_mask):
                    reprojected_data[nan_mask] = 0
                    print(f"      Replaced {np.sum(nan_mask)} NaN pixels with zeros")
                
                # Create new FITS file with reference WCS
                new_header = reference_wcs.to_header()
                
                # Preserve important original header information
                preserve_keys = ['OBJECT', 'FILTER', 'EXPTIME', 'DATE-OBS', 'TELESCOP', 'INSTRUME', 
                                'BUNIT', 'COMMENT', 'HISTORY']
                for key in preserve_keys:
                    if key in header:
                        if key in ['COMMENT', 'HISTORY']:
                            # Handle multiple entries
                            for value in header[key]:
                                new_header[key] = value
                        else:
                            new_header[key] = header[key]
                
                # Add info about reprojection
                new_header['HISTORY'] = f'Reprojected to common WCS grid by create_fitsmap.py'
                new_header['HISTORY'] = f'Original file: {os.path.basename(input_file)}'
                
                # Create primary HDU with reprojected data
                primary_hdu = fits.PrimaryHDU(data=reprojected_data, header=new_header)
                hdul_new = fits.HDUList([primary_hdu])
                
                # Save aligned file
                hdul_new.writeto(output_file, overwrite=True)
                print(f"      Successfully saved aligned version to {os.path.basename(output_file)}")
                return True
                
            except Exception as e:
                print(f"      Error during reprojection for {input_file}: {e}")
                print(f"      Attempting to copy original file as fallback...")
                # Copy the original file if reprojection fails
                with open(input_file, 'rb') as src, open(output_file, 'wb') as dst:
                    dst.write(src.read())
                print(f"      Copied original file (WCS offset may remain)")
                return False
                
    except Exception as e:
        print(f"      Error reading {input_file}: {e}")
        return False

def create_fitsmap(input_dir, output_dir, reference_file=None, serve=False, colormap="gray", align_images=True, add_coord_search=True, debug=False):
    """
    Create FITSMap tiles from preprocessed FITS files and RGB images.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing preprocessed FITS files and RGB images
    output_dir : str
        Directory where FITSMap will be generated
    reference_file : str, optional
        Specific file to use as WCS reference
    serve : bool
        Whether to start the web server after generation
    colormap : str
        Matplotlib colormap to use for rendering FITS files
    align_images : bool
        Whether to align FITS images to common WCS grid
    add_coord_search : bool
        Whether to add coordinate search functionality to the viewer
    """
    # Find all FITS and image files
    fits_files, image_files = find_fits_and_image_files(input_dir)
    
    if not fits_files and not image_files:
        print(f"No FITS or image files found in {input_dir}")
        return
    
    print(f"Found {len(fits_files)} FITS files:")
    for f in fits_files:
        print(f"  - {os.path.basename(f)}")
    
    if image_files:
        print(f"Found {len(image_files)} image files:")
        for f in image_files:
            print(f"  - {os.path.basename(f)}")
    
    # Check WCS alignment if we have multiple FITS files
    if len(fits_files) > 1:
        check_wcs_alignment(fits_files)
    
    # Handle FITS alignment if requested and multiple FITS files exist
    working_dir = input_dir
    temp_dir = None
    
    if align_images and len(fits_files) > 1:
        print(f"\nAligning FITS images to common WCS grid...")
        
        # Find reference WCS
        if reference_file is None:
            ref_file, ref_wcs, ref_shape = find_reference_wcs(fits_files)
        else:
            ref_wcs, ref_shape, _ = get_wcs_and_shape(reference_file)
            ref_file = reference_file
        
        if ref_wcs is None:
            print("Warning: Could not determine reference WCS, skipping alignment")
        else:
            # Create temporary directory for aligned files
            temp_dir = tempfile.mkdtemp(prefix='fitsmap_aligned_')
            print(f"Creating aligned files in temporary directory: {temp_dir}")
            
            # Copy image files to temp directory
            for img_file in image_files:
                temp_img = os.path.join(temp_dir, os.path.basename(img_file))
                with open(img_file, 'rb') as src, open(temp_img, 'wb') as dst:
                    dst.write(src.read())
                print(f"    Copied {os.path.basename(img_file)} to temp directory")
            
            # Check all FITS files for WCS offsets and reproject them
            aligned_fits = []
            reprojection_needed = False
            
            print(f"\nProcessing {len(fits_files)} FITS files:")
            for fits_file in fits_files:
                temp_fits = os.path.join(temp_dir, os.path.basename(fits_file))
                
                # Always check for WCS differences, even for the reference file
                file_wcs, file_shape, _ = get_wcs_and_shape(fits_file)
                
                if file_wcs is not None and ref_wcs is not None:
                    # Compare WCS centers
                    file_center = file_wcs.wcs.crval if hasattr(file_wcs.wcs, 'crval') else None
                    ref_center = ref_wcs.wcs.crval if hasattr(ref_wcs.wcs, 'crval') else None
                    
                    if file_center is not None and ref_center is not None:
                        ra_offset = file_center[0] - ref_center[0]
                        dec_offset = file_center[1] - ref_center[1]
                        total_offset_arcsec = np.sqrt(ra_offset**2 + dec_offset**2) * 3600
                        
                        print(f"  {os.path.basename(fits_file)}:")
                        print(f"    Center: RA={file_center[0]:.6f}°, Dec={file_center[1]:.6f}°")
                        print(f"    Offset from reference: {total_offset_arcsec:.2f} arcsec")
                        
                        # If this is the reference file, just copy it
                        if fits_file == ref_file:
                            print(f"    This is the reference file - copying as-is")
                            with open(fits_file, 'rb') as src, open(temp_fits, 'wb') as dst:
                                dst.write(src.read())
                        else:
                            # Always reproject non-reference files to ensure alignment
                            print(f"    Reprojecting to match reference WCS...")
                            success = reproject_fits_file(fits_file, ref_wcs, ref_shape, temp_fits)
                            if success:
                                reprojection_needed = True
                            else:
                                print(f"    Warning: Reprojection failed for {fits_file}")
                    else:
                        print(f"  {os.path.basename(fits_file)}: Could not determine WCS center")
                        # Copy without reprojection
                        with open(fits_file, 'rb') as src, open(temp_fits, 'wb') as dst:
                            dst.write(src.read())
                else:
                    print(f"  {os.path.basename(fits_file)}: No valid WCS found")
                    # Copy without reprojection
                    with open(fits_file, 'rb') as src, open(temp_fits, 'wb') as dst:
                        dst.write(src.read())
                
                aligned_fits.append(temp_fits)
            
            # Use temp directory as working directory
            working_dir = temp_dir
            reference_file = os.path.join(temp_dir, os.path.basename(ref_file))
            
            if reprojection_needed:
                print(f"\nSuccessfully aligned files to {os.path.basename(ref_file)} WCS grid")
            else:
                print(f"\nUsing {os.path.basename(ref_file)} as reference (no reprojection needed)")
    elif align_images and len(fits_files) == 1:
        print(f"\nOnly one FITS file found - no alignment needed")
    elif not align_images:
        print(f"\nWCS alignment disabled (--no-align specified)")
    else:
        print(f"\nNo FITS files found for alignment")
    
    # Select WCS reference file (must be a FITS file)
    if reference_file is None:
        if fits_files:
            reference_file = select_wcs_reference(fits_files)
        else:
            print("Error: No FITS files found for WCS reference")
            print("PNG files do not contain WCS information - need at least one FITS file")
            if temp_dir:
                import shutil
                shutil.rmtree(temp_dir)
            return
    
    if reference_file is None:
        print("Error: No suitable WCS reference file found")
        if temp_dir:
            import shutil
            shutil.rmtree(temp_dir)
        return
    
    print(f"Using WCS reference: {os.path.basename(reference_file)}")
    
    # Set colormap for FITSMap (applies to FITS files only)
    convert.MPL_CMAP = colormap
    print(f"Using colormap for FITS files: {colormap}")
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Normalization settings
    norm_kwargs = {
        "stretch": "linear",
        "min_percent": 0,
        "max_percent": 100,
    }
    
    print(f"Generating FITSMap in: {output_dir}")
    print("FITS normalization: Linear stretch, no additional clipping")
    if image_files:
        print("PNG/Image files: Displayed as-is")
    if align_images and len(fits_files) > 1:
        print("FITS files: Aligned to common WCS grid")
    print("(preserving preprocessed normalization)")
    
    try:
        # Generate the FITSMap using aligned files
        convert.dir_to_map(
            working_dir,
            out_dir=output_dir,
            cat_wcs_fits_file=reference_file,
            norm_kwargs=norm_kwargs
        )
        
        print("FITSMap successfully generated!")
        print(f"Map location: {output_dir}")
        print(f"Main file: {os.path.join(output_dir, 'index.html')}")
        
        if fits_files:
            alignment_status = " (WCS-aligned)" if align_images and len(fits_files) > 1 else ""
            print(f"FITS layers: {len(fits_files)} files{alignment_status}")
        if image_files:
            print(f"Image layers: {len(image_files)} files (including your RGB image)")
        
        # Add coordinate search functionality
        if add_coord_search:
            add_coordinate_search(output_dir, reference_file)
        
        if serve:
            print("Starting web server...")
            serve_fitsmap(output_dir)
        else:
            print("To view the map, run:")
            print(f"  cd {output_dir}")
            print(f"  fitsmap serve")
            print("Or use this script with --serve flag")
            
    except Exception as e:
        print(f"Error generating FITSMap: {e}")
        raise
    finally:
        # Clean up temporary directory
        if temp_dir:
            print("Cleaning up temporary aligned files...")
            import shutil
            shutil.rmtree(temp_dir)

def serve_fitsmap(map_dir):
    """
    Start the FITSMap web server.
    
    Parameters:
    -----------
    map_dir : str
        Directory containing the generated FITSMap
    """
    try:
        # Change to map directory
        original_dir = os.getcwd()
        os.chdir(map_dir)
        
        print(f"Starting FITSMap server in {map_dir}")
        print("Press Ctrl+C to stop the server")
        
        # Start the server
        subprocess.run(["fitsmap", "serve"])
        
    except KeyboardInterrupt:
        print("Stopping server...")
    except FileNotFoundError:
        print("Error: 'fitsmap' command not found.")
        print("Make sure FITSMap is installed and in your PATH.")
        print("You can manually serve by running 'fitsmap serve' in the map directory.")
    finally:
        # Return to original directory
        os.chdir(original_dir)

def check_fitsmap_installation():
    """Check if FITSMap is properly installed."""
    try:
        import fitsmap
        # Try to get version, but don't fail if it doesn't exist
        try:
            version = fitsmap.__version__
            print(f"FITSMap version: {version}")
        except AttributeError:
            print("FITSMap found (version unknown)")
        return True
    except ImportError:
        print("Error: FITSMap not found.")
        print("Install with: pip install fitsmap")
        return False

def main():
    """Main function with command line interface."""
    parser = argparse.ArgumentParser(
        description="Generate FITSMap tiles from preprocessed JWST FITS files and RGB images",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate FITSMap with WCS-aligned layers and coordinate search (recommended)
  python create_fitsmap.py /path/to/data /path/to/output/map

  # Generate and immediately serve
  python create_fitsmap.py /path/to/data /path/to/output/map --serve

  # Skip alignment for faster processing (may have layer offsets)
  python create_fitsmap.py /path/to/data /path/to/output/map --no-align

  # Use specific WCS reference file (if automatic selection is wrong)
  python create_fitsmap.py input_dir output_dir --reference PSZ2G147_f150w2_i2d.fits

  # Debug WCS alignment issues
  python create_fitsmap.py input_dir output_dir --debug

  # Use different colormap for FITS files
  python create_fitsmap.py input_dir output_dir --colormap viridis

Troubleshooting WCS Offsets:
  If you see layer offsets in the FITSMap viewer:
  1. Run with --debug to see detailed WCS information
  2. Try specifying a different reference file with --reference
  3. Check that your FITS files have valid WCS headers
  4. The script will show offset calculations before alignment

Input Directory Structure:
  /path/to/data/
    ├── PSZ2G147_f322w2_i2d.fits  # FITS files (for individual bands)
    ├── PSZ2G147_f150w2_i2d.fits  
    ├── rgb_image.png             # Your RGB composite image
    └── catalog.cat               # Optional: source catalog

The script will automatically:
  - Detect and report WCS offsets between FITS files
  - Align all FITS files to a common WCS grid (eliminates offsets)
  - Create individual grayscale layers for each FITS file
  - Include your RGB image as a color layer  
  - Add coordinate search box (bottom-left corner)
  - Use proper WCS coordinates for all layers
  - Enable layer switching/blending in the browser interface

Features:
  WCS Alignment: Aligns all FITS files to finest-resolution WCS grid
  Coordinate Search: Enter RA/Dec coordinates to center the map
    - Supports both decimal degrees (150.5, -20.3) 
    - And sexagesimal format (10:02:00, -20:18:00)
    - Box appears in bottom-left corner of the viewer
        """
    )
    
    parser.add_argument('input_dir', help='Directory containing FITS files and RGB images')
    parser.add_argument('output_dir', help='Output directory for FITSMap')
    parser.add_argument('--reference', help='Specific FITS file to use as WCS reference (full path or filename in input_dir)')
    parser.add_argument('--serve', action='store_true', help='Start web server after generation')
    parser.add_argument('--colormap', default='gray', help='Matplotlib colormap for FITS files (default: gray)')
    parser.add_argument('--no-align', action='store_true', help='Skip WCS alignment of FITS files (faster but may have offsets)')
    parser.add_argument('--no-coord-search', action='store_true', help='Skip adding coordinate search functionality')
    parser.add_argument('--debug', action='store_true', help='Enable detailed debugging output for WCS alignment')
    
    args = parser.parse_args()
    
    # Check if FITSMap is installed
    if not check_fitsmap_installation():
        sys.exit(1)
    
    # Check if input directory exists
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' does not exist")
        sys.exit(1)
    
    if not os.path.isdir(args.input_dir):
        print(f"Error: '{args.input_dir}' is not a directory")
        sys.exit(1)
    
    # Validate reference file if provided
    if args.reference:
        if os.path.exists(args.reference):
            # Full path provided
            reference_file_path = args.reference
        else:
            # Check if it's a filename within input_dir
            reference_file_path = os.path.join(args.input_dir, args.reference)
            if not os.path.exists(reference_file_path):
                print(f"Error: Reference file '{args.reference}' not found")
                print(f"Checked: {args.reference}")
                print(f"Checked: {reference_file_path}")
                sys.exit(1)
        args.reference = reference_file_path
    
    # Generate FITSMap
    create_fitsmap(
        args.input_dir,
        args.output_dir,
        reference_file=args.reference,
        serve=args.serve,
        colormap=args.colormap,
        align_images=not args.no_align,
        add_coord_search=not args.no_coord_search,
        debug=args.debug
    )

if __name__ == "__main__":
    main()