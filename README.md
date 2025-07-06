# Interpolation Project

This repository contains a Python script (`interpolate_and_plot.py`) for interpolating velocity fields and generating slice plots. It demonstrates the use of `numpy`, `scipy`, `matplotlib`, and `vtk` for data generation, interpolation, visualization, and VTK export.

## Setup

To set up the project, it is recommended to use a virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate
pip install numpy scipy matplotlib vtk
```

## Usage

To run the script and generate the `velocity_data.vtu` file and `velocity_difference_slice.png` plot:

```bash
source venv/bin/activate
python interpolate_and_plot.py
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.