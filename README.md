# Interpolation Project

This repository contains a Python script (`interpolate_and_plot.py`) for interpolating velocity fields and generating slice plots. It demonstrates the use of `numpy`, `scipy`, `matplotlib`, and `vtk` for data generation, interpolation, visualization, and VTK export. 
This repo currently uses a demo with random arrays. The data that is meant to be used with it is gnerated by `MUSIC` version 2 (https://github.com/cosmo-sims/MUSIC2)

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

## Demo

![Velocity Difference Slice Plot](https://github.com/user-attachments/assets/f3fb52b2-32b0-4eba-9afd-006dc90b4683)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
