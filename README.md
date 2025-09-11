# slice-data-viz

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

Data Visualizer for RGB images of JWST-SLICE Galaxy Clusters

## Overview

`slice-data-viz` is a specialized data visualization tool designed for analyzing and visualizing RGB images of galaxy clusters observed by the James Webb Space Telescope (JWST) as part of the SLICE (Surveys of Lensing Clusters and Environments) program. This tool provides researchers and astronomers with powerful capabilities to process, analyze, and create compelling visualizations of distant galaxy clusters.

## Features

- **JWST Image Processing**: Specialized tools for handling JWST RGB image data
- **Galaxy Cluster Analysis**: Advanced algorithms for analyzing galaxy cluster structures
- **Interactive Visualizations**: Create dynamic and interactive plots for data exploration
- **RGB Image Enhancement**: Tools for enhancing and processing RGB astronomical images
- **Data Export**: Multiple export formats for sharing and publication
- **Customizable Plotting**: Flexible visualization options tailored for astronomical data

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Quick Install

```bash
git clone https://github.com/gkhullar/slice-data-viz.git
cd slice-data-viz
pip install -r requirements.txt
```

### Development Installation

```bash
git clone https://github.com/gkhullar/slice-data-viz.git
cd slice-data-viz
pip install -e .
```

## Usage

### Basic Usage

```python
import slice_data_viz as sdv

# Load JWST RGB image data
data = sdv.load_jwst_image('path/to/jwst_image.fits')

# Create visualization
viz = sdv.GalaxyClusterVisualizer(data)
viz.plot_rgb_image()
viz.show()
```

### Advanced Features

```python
# Analyze galaxy cluster properties
cluster = sdv.ClusterAnalysis(data)
cluster.identify_members()
cluster.calculate_mass_distribution()

# Create publication-ready plots
plotter = sdv.PublicationPlotter(cluster)
plotter.create_figure(style='publication')
plotter.save('galaxy_cluster_analysis.png', dpi=300)
```

## Data Requirements

This tool is designed to work with:
- JWST RGB image data (FITS format)
- Galaxy cluster photometric catalogs
- Redshift data for cluster members
- Multi-wavelength astronomical imaging data

## Contributing

We welcome contributions to `slice-data-viz`! Please see our contributing guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Development Setup

```bash
# Clone your fork
git clone https://github.com/yourusername/slice-data-viz.git
cd slice-data-viz

# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/

# Run linting
flake8 slice_data_viz/
black slice_data_viz/
```

## Documentation

Detailed documentation is available at [documentation link] (coming soon).

For API reference and examples, see the `docs/` directory.

## Citation

If you use `slice-data-viz` in your research, please cite:

```bibtex
@software{slice_data_viz,
  author = {Khullar, Gourav},
  title = {slice-data-viz: Data Visualizer for RGB images of JWST-SLICE Galaxy Clusters},
  url = {https://github.com/gkhullar/slice-data-viz},
  year = {2025}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- James Webb Space Telescope team
- SLICE collaboration
- Astronomical data visualization community

## Support

For questions, issues, or feature requests, please:
- Open an issue on GitHub
- Contact the maintainer: [Gourav Khullar](https://github.com/gkhullar)

## Related Projects

- [JWST Pipeline](https://github.com/spacetelescope/jwst)
- [Astropy](https://www.astropy.org/)
- [Matplotlib Astronomy](https://matplotlib.org/stable/gallery/index.html#astronomy)

---

**Note**: This project is under active development. Features and API may change as the project evolves.
