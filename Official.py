# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 16:00:56 2024

@author: mike_
"""
from lightkurve import TessTargetPixelFile
import numpy as np

def process_tpf(file_path, period_range=(1, 20), window_length=901, frequency_factor=500, max_planets=2, plot_xlim=(-1, 1)):
    """
    Function to process a TESS Target Pixel File, search for multiple transit signals, and plot results.
    
    Parameters:
    - file_path: Path to the target pixel file.
    - period_range: Tuple (min_period, max_period) to define the range of periods for the BLS search.
    - window_length: Length of the window used to flatten the light curve.
    - frequency_factor: Factor used to calculate the number of frequencies for BLS search.
    - max_planets: Maximum number of planets to detect.
    - plot_xlim: Range for x-axis when plotting folded light curves.
    """
    # Load and plot the target pixel file
    tpf = TessTargetPixelFile(file_path)
    tpf.plot()
    
    # Convert to light curve, flatten it and remove outliers
    lc = tpf.to_lightcurve().flatten(window_length=window_length).remove_outliers()
    lc.plot()
    
    # Create an array of periods to search for the planet's transit signal
    period = np.linspace(period_range[0], period_range[1], 10000)
    
    # Initialize list to store detected planet parameters
    detected_planets = []
    
    # Iteratively search for planets
    for i in range(max_planets):
        # Perform a BLS search to find the best-fitting transit model
        bls = lc.to_periodogram(method='bls', period=period, frequency_factor=frequency_factor)
        bls.plot()
        
        # Get the best-fit transit parameters (period, transit time, duration)
        planet_period = bls.period_at_max_power
        planet_t0 = bls.transit_time_at_max_power
        planet_dur = bls.duration_at_max_power
        
        # Store the detected planet parameters
        detected_planets.append((planet_period, planet_t0, planet_dur))
        print(f"Detected Planet {i + 1} Period: {planet_period}")
        
        # Create the BLS transit model using the detected parameters
        planet_model = bls.get_transit_model(period=planet_period, transit_time=planet_t0, duration=planet_dur)
        
        # Fold the light curve on the detected period and plot it with the BLS model
        ax = lc.fold(planet_period, planet_t0).scatter(s=0.5)
        planet_model.fold(planet_period, planet_t0).plot(ax=ax, c='r', lw=2, label=f'Transit Model {i + 1}')
        ax.set_xlim(plot_xlim)
        ax.legend()
        
        # Subtract the detected planet's signal from the light curve (prewhitening)
        lc = lc - planet_model
        

    
    return detected_planets


# Example 1:
process_tpf("tess2020212050318-s0028-0000000266980320-0190-a_fast-tp.fits", max_planets=1)

# Example 2:
process_tpf("tess2018234235059-s0002-0000000270622440-0121-s_tp.fits")

# Example 3: Specify a different period range and x-axis limits
process_tpf("tess2019279210107-s0017-0000000445258206-0161-s_tp.fits", period_range=(1.5, 20), plot_xlim=(-0.5, 0.5),max_planets=1)

# Example 4:
process_tpf("tess2018206045859-s0001-0000000263003176-0120-s_tp.fits",max_planets=1)

# Example 5:
process_tpf("tess2018206045859-s0001-0000000052368076-0120-s_tp.fits",max_planets=1)

