# SLAM Mirror Stellarator Campaign

This directory contains reusable mirror-campaign simulation tools, analysis scripts, and results from the WarpX hackathon beam anisotropy project.

## Contents

- **inputs/**: Example input files for SLAM mirror magnetic confinement simulations
- **scripts/**: Python generators, analyzers, and visualization tools for mirror campaigns
- **summaries/**: Compact CSV and JSON summaries of campaign results, comparisons, and validation metrics
- **figures/**: Selected PNG visualizations and preview images from campaign analysis

## Scripts

- `generate_inputs.py` - Generate simulation input files from parameter templates
- `generate_campaign.py` - Create a suite of parameter-varied simulations
- `analyze_case.py` - Analyze a single simulation result
- `analyze_campaign.py` - Batch analyze campaign results
- `plot_campaign.py` - Create summary plots of campaign data
- `create_visualizations.py` - Generate detailed analysis figures
- `animate_campaign_3d.py` - Create 3D trajectory animations
- `animate_campaign_velocity_space.py` - Create velocity-space evolution animations
- `assess_horizon.py` - Evaluate confinement and loss cone metrics
- `compare_angle_pilot.py` - Compare pilot and full-angle campaigns
- `summarize_angle_pilot.py` - Produce summary statistics
- `summarize_validation.py` - Validate simulation against reference results

## Usage

These tools are designed to integrate with WarpX simulations and support parameter studies, sensitivity analysis, and visualization of charged particle confinement in magnetic mirror geometries.

## References

Generated during the Plasma Hackathon 2024 hackathon-beam-anisotropy project.
