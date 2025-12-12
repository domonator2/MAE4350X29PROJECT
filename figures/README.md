# Figures

This directory contains all generated plots and figures from the analysis.

## Contents

All 43 PNG figures from the trade study analysis, including:

### Design Trade Study Plots
- Carpet plots (TOGW, OEW, fuel weight, structural weight, systems weight)
- Wing loading analysis
- T/W ratio analysis
- Solution space visualizations

### Aerodynamic Analysis
- L/D polars
- Drag breakdown
- Mach sweep analysis

### Convergence Analysis
- TOGW vs sweep angle
- TOGW vs slenderness (tau)
- Fuel fraction vs duration
- Geometry constraints

### Comparison Plots
- Raymer vs Roskam weight methods
- Elliptical fuselage trade studies
- Sensitivity analyses

## Usage

Figures are referenced in:
- Final project report
- Presentations
- Technical documentation
- Trade study analysis

## Organization

Figures follow the naming convention:
- `carpet_plot_*.png` - Carpet plot visualizations
- `*_vs_*.png` - Comparison/trade plots
- `*_analysis.png` - Analysis results
- `*_comparison.png` - Method comparisons

## Generation

Most figures are generated automatically by scripts in:
- `Synthesis/finalized_convergence/execute_design_trade_study.m`
- `test_debug/plotting/` scripts
