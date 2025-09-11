# slice-data-viz

**Data Visualizer for RGB images of JWST-SLICE Galaxy Clusters**
!([http://url/to/img.png](https://github.com/gkhullar/slice-data-viz/blob/main/Logo%20slice.jpeg))

This repository provides tools for processing James Webb Space Telescope (JWST) observations of galaxy clusters from the SLICE survey and creating beautiful RGB composite images. The pipeline handles FITS file preprocessing, coordinate alignment, and advanced RGB image creation with astronomy-specific stretch algorithms.

## Features

- **FITS File Processing**: Handles JWST multi-extension FITS files with automatic SCI extension detection
- **DS9-Style Normalization**: Applies asinh stretching with configurable contrast and bias parameters
- **Color Inversion**: Optional astronomical convention color inversion (black background, white stars)
- **Multi-Filter RGB Composition**: Creates RGB images from different JWST filter observations
- **Batch Processing**: Processes entire directory structures of galaxy cluster observations
- **WCS-Aware Alignment**: Automatic coordinate system alignment using World Coordinate Systems
- **Advanced Stretch Algorithms**: TRILOGY-style RGB creation with per-channel parameter control
- **Flexible Filter Mapping**: Case-insensitive filter identification for RGB channel assignment

## Installation

### Prerequisites

This tool requires Python 3.7+ and several astronomical and image processing libraries:

```bash
pip install astropy reproject matplotlib numpy Pillow
```

### Clone the Repository

```bash
git clone https://github.com/gkhullar/slice-data-viz.git
cd slice-data-viz
```

## Directory Structure

```
slice-data-viz/
├── README.md                          # This file
├── LICENSE                            # License file
├── combined_fits_rgb_processor.py     # Main processing script
├── PSZ2G147_RGB.png                  # Example RGB output image
├── data/                             # Input JWST FITS files
│   └── readme.md                     # Data directory instructions
├── processed_fits/                   # Processed/normalized FITS outputs
│   └── readme.md
├── output_maps/                      # Additional output directory
│   └── readme.md
└── diagnostics_and_reserve/         # Diagnostic outputs and backup files
    └── readme.md
```

## Quick Start

### 1. Prepare Your Data

Place your JWST FITS files in subdirectories under `data/`, organized by galaxy cluster:

```
data/
├── cluster_name_1/
│   ├── observation_F150W2.fits
│   ├── observation_F322W2.fits
│   └── observation_F444W.fits
└── cluster_name_2/
    ├── observation_F150W2.fits
    └── observation_F322W2.fits
```

### 2. Basic Usage

Process all clusters with default settings:

```bash
python combined_fits_rgb_processor.py --input-dir ./data --output-dir ./processed_fits_rgb
```

### 3. Custom Filter Mapping

Specify which filters to use for RGB channels:

```bash
python combined_fits_rgb_processor.py \
    --input-dir ./data \
    --output-dir ./processed_fits_rgb \
    -r f322w2 -g f322w2 -b f150w2
```

## Usage Examples

### Default Processing
```bash
# Process with default filter mapping (f322w2 for red/green, f150w2 for blue)
python combined_fits_rgb_processor.py --input-dir ./data --output-dir ./output
```

### Custom Filter Combinations
```bash
# Use different filter combinations (case-insensitive)
python combined_fits_rgb_processor.py \
    --input-dir ./data --output-dir ./output \
    -r F444W -g F277W -b F150W
```

### Advanced RGB Parameter Tuning
```bash
# Fine-tune RGB creation with per-channel parameters
python combined_fits_rgb_processor.py \
    --input-dir ./data --output-dir ./output \
    -r f322w2 -g f322w2 -b f150w2 \
    --b-asinh-softening 0.1 \
    --b-brightness-boost 1.5 \
    --b-noiselum 0.3
```

### Conservative Preprocessing
```bash
# Use conservative preprocessing parameters
python combined_fits_rgb_processor.py \
    --input-dir ./data --output-dir ./output \
    --contrast 1.5 --bias 0.5 --max-percent 95
```

### High Contrast RGB
```bash
# Create high-contrast images
python combined_fits_rgb_processor.py \
    --input-dir ./data --output-dir ./output \
    --asinh-softening 0.02 --brightness-boost 5.0
```

## Command Line Options

### Input/Output
- `--input-dir`: Input directory containing galaxy subdirectories
- `--output-dir`: Output directory for processed files and RGB images

### RGB Channel Mapping (case-insensitive)
- `-r, --red`: Filter code for red channel (default: f322w2)
- `-g, --green`: Filter code for green channel (default: f322w2)  
- `-b, --blue`: Filter code for blue channel (default: f150w2)

### FITS Preprocessing Parameters
- `--contrast`: DS9 contrast parameter (default: 1.02)
- `--bias`: DS9 bias parameter (default: 0.43)
- `--max-percent`: Maximum percentile for data clipping (default: 99)

### RGB Creation Parameters
- `--satpercent`: Percentage of pixels to saturate (default: 0.001)
- `--colorsatfac`: Color saturation factor (default: 1)
- `--noiselum`: Noise luminosity threshold (default: 0.15)
- `--asinh-softening`: Default asinh softening parameter (default: 0.05)
- `--brightness-boost`: Default brightness multiplier (default: 3.0)

### Per-Channel Parameters
Individual control for each RGB channel:
- `--r-asinh-softening`, `--g-asinh-softening`, `--b-asinh-softening`
- `--r-brightness-boost`, `--g-brightness-boost`, `--b-brightness-boost`
- `--r-noiselum`, `--g-noiselum`, `--b-noiselum`
- `--r-scale`, `--g-scale`, `--b-scale`

## Technical Details

### FITS Processing Pipeline

1. **Extension Detection**: Automatically finds SCI extensions in JWST multi-extension FITS files
2. **DS9 Normalization**: Applies asinh stretching with configurable parameters
3. **Color Inversion**: Converts astronomical data to conventional RGB display format
4. **Header Updates**: Preserves metadata and adds processing history

### RGB Creation Process

1. **WCS Alignment**: Uses World Coordinate System information for accurate alignment
2. **Reprojection**: Aligns images to common coordinate grid using `reproject` library
3. **TRILOGY Algorithm**: Implements advanced RGB stretch algorithms from TRILOGY software
4. **Per-Channel Control**: Allows independent parameter tuning for each color channel

### Supported JWST Filters

The tool works with any JWST filter observations. Common filter combinations:
- **Near-Infrared**: F150W, F200W, F277W, F356W, F444W
- **Wide Filters**: F150W2, F322W2 (commonly used in SLICE survey)

## Data Requirements

### Input Format
- JWST FITS files with SCI extensions
- Files organized in subdirectories by galaxy cluster
- Filter information embedded in filenames (case-insensitive matching)

### Naming Convention
Files should include filter names in their filenames for automatic detection:
- `cluster_observation_F150W2.fits`
- `PSZ2G147_f322w2_processed.fits`
- Any filename containing the filter code (case-insensitive)

## Output

For each processed galaxy cluster, the pipeline generates:

1. **Processed FITS Files**: Normalized, single-extension FITS files in `output_dir/cluster_name/`
2. **RGB Images**: High-quality PNG images named `cluster_name_RGB.png`
3. **Processing Metadata**: FITS headers contain processing history and parameters used

## Example Output

See `PSZ2G147_RGB.png` for an example of the RGB composite image quality achievable with this pipeline.

## Contributing

When contributing to this project:

1. Ensure all dependencies are properly documented
2. Test with different JWST filter combinations
3. Maintain compatibility with astronomical FITS standards
4. Document any new command-line parameters

## Dependencies

- **astropy**: FITS file handling and astronomical coordinate systems
- **reproject**: Image reprojection and alignment
- **matplotlib**: Color mapping and image processing utilities  
- **numpy**: Numerical computations and array operations
- **Pillow (PIL)**: Final image creation and file output

## License

See `LICENSE` file for license information.

## Related Projects

This tool is designed for the SLICE (Survey of Lensing clusters In CosmoS-simulation) galaxy cluster survey using JWST observations.
