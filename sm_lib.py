import matplotlib.pyplot as plt
import numpy as np
import uproot
import sifi_cm.root_aux as raux
from sifi_cm import distal
from sifi_cm.data_fit import smooth, fit_1d, Gaussian_1D, normalize
from pathlib import Path
import mplhep
from scipy.stats import pearsonr
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import re
from scipy.stats import pearsonr
from numpy.linalg import norm

plt.style.use(mplhep.style.LHCb2)

def reco_mlem_from_sim(dataPath, histName, matrFile, nIter, xRange, 
         rangesPSTAR=None, compare_pstar=True, only_ranges=False, verbose=False):
    """
    Reconstructs histograms from ROOT files using MLEM,
    detects peaks, and (optionally) compares to PSTAR reference ranges.

    Notes
    -----
    - ROOT files inside `dataPath` must follow one of the formats:
        1. "run00{runid}_*.root"    → key = run ID (int)
        2. "{energy}_*.root" or "{energy}dot{fraction}_*.root"       → key = energy in MeV (str)

      If the filename contains '_bg', the file is skipped.

    Parameters
    ----------
    dataPath : str or Path
        Directory with ROOT files.
    histName : str
        Histogram name inside ROOT files.
    matrFile : ndarray
        System matrix used for MLEM.
    nIter : int
        Number of MLEM iterations.
    xRange : array-like
        Range values (mm) corresponding to histogram bins.
    rangesPSTAR : dict, optional
        PSTAR reference ranges keyed by runID or energy.
    compare_pstar : bool, default=True
        Whether to perform PSTAR comparison.
    only_ranges : bool, default=False
        If True, skip plotting and only return ranges.
    verbose : bool, default=False
        If True, prints detailed info.

    Returns
    -------
    result : dict
        {
          "key_type": "runid" or "energy",
          "ranges": {key: detected_range_mm, ...},
          "foundRange": list of detected peaks [mm],
          "fig": matplotlib.figure.Figure or None
        }
    """
    dataPath = Path(dataPath)

    rangeList = {}   # detected peak positions keyed by run ID / energy
    foundRange = []  # list of all detected peak positions
    key_type = None  # "runid" or "energy"

    fig, axs = (None, None)
    if not only_ranges:
        fig, axs = plt.subplots(2, 1, figsize=(20, 15), 
                                 gridspec_kw={'height_ratios': [2, 1]})
        cmap = plt.get_cmap("tab20")

    # --- Process each ROOT file in directory ---
    profiles = []
    for i, filePath in enumerate(sorted(dataPath.glob("*.root"))):

        parts = filePath.stem.split('_')

        # Skip background files
        if any("bg" in p.lower() for p in parts):
            if verbose:
                print(f"Skipping background file: {filePath.name}")
            continue

        # Determine key type
        if parts[0].startswith("run"):
            if key_type is None:
                key_type = "runid"
            key = int(parts[0].lstrip("run0"))
            if verbose:
                print(f"Processing runID: {key}")
        else:
            if key_type is None:
                key_type = "energy"
            key = str(parts[0].rstrip("_"))
            if verbose:
                print(f"Processing energy: {key} MeV")

        if not filePath.is_file():
            continue

        try:
            # --- Load histogram ---
            simHist = raux.get_histo(path=filePath, histo_names=histName)
            image = simHist.vals.T[::-1]

            # --- Reconstruct using MLEM ---
            
            matr = matrFile.T[::-1].T
            
               
            sens = np.sum(matr)
            simReco = raux.reco_mlem(matr=matr, image=image.flatten(), niter=nIter, S=sens)

            # --- Apply Gaussian smoothing ---
            simReco_smoothed = smooth(simReco, scale=3, filter='gaussian')

            # --- Peak detection ---
            peaks, _ = find_peaks(simReco_smoothed, height=0.1, prominence=0.05)

            if peaks.size > 0:
                peak_idx = peaks[0]
                if verbose:
                    print(f"{key_type} {key}: peak at {xRange[::-1][peak_idx]:.2f}")
            else:
                peak_idx = np.argmax(simReco_smoothed)
                if verbose:
                    print(f"No clear peak for {key_type} {key}, using max at {xRange[::-1][peak_idx]:.2f}")


            # --- Save detected peak position ---
            rangeList[key] = xRange[::-1][peak_idx]
            foundRange.append(xRange[::-1][peak_idx])

                    # Convert energy string to float
            if key_type == "energy":
                # Remove non-numeric parts and handle "21dot1" -> 21.1
                energy_value = float(str(key).replace("dot", "."))
            else:
                energy_value = key

            profiles.append((energy_value, simReco_smoothed, peak_idx))

        except uproot.KeyInFileError as e:
            print(f"⚠️ Skipping {filePath} - Missing key: {e}")

    if only_ranges and not compare_pstar:
        return {
            "key_type": key_type,
            "ranges": rangeList,
            "foundRange": foundRange,
            
        }
    if not only_ranges:
        # Sort profiles by energy ascending
        profiles.sort(key=lambda x: x[0])

        # Plot
        cmap = plt.get_cmap("tab20")
        for i, (energy_value, profile, peak_idx) in enumerate(profiles):
            color = cmap(i % 20)
            label = f"{energy_value:.1f} MeV" if key_type == "energy" else str(energy_value)
            axs[0].plot(xRange[::-1], profile.reshape(100, 1), label=label, color=color)
            axs[0].axvline(xRange[::-1][peak_idx], color=color, linestyle='--', linewidth=2)

        axs[0].legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)



    # --- Optional PSTAR comparison ---
    if compare_pstar and rangesPSTAR is not None:
        expected_ranges = []
        observed_ranges = []

        for run_id in rangeList:
            if run_id in rangesPSTAR:
                if verbose:
                    print(f'PSTAR = {rangesPSTAR[run_id]}, Peak = {rangeList[run_id]}')
                expected_ranges.append(rangesPSTAR[run_id])
                observed_ranges.append(rangeList[run_id])
        

        coef, cov = np.polyfit(expected_ranges, observed_ranges, 1, cov=True)
        poly_fn = np.poly1d(coef)
        slope_err = np.sqrt(cov[0, 0])
        expected = np.array(expected_ranges)
        observed = np.array(observed_ranges)
        y_pred = poly_fn(expected)

        rmse = np.sqrt(np.mean((observed - y_pred) ** 2))
        corr_coef, p_value = pearsonr(expected, observed)

        if not only_ranges:
            axs[1].scatter(expected_ranges, observed_ranges, label='Data')
            axs[1].plot(expected_ranges, poly_fn(expected_ranges),
                    label=f'slope = {coef[0]:.3f} ± {slope_err:.3f}', color='red')

            axs[1].legend()
            axs[1].set_ylabel('Observed peak [mm]')
            axs[1].set_xlabel('PSTAR expected [mm]')
            axs[1].text(0.05, 0.8, f"RMSE = {rmse:.4f}\nCorr. coeff = {corr_coef:.4f}",
                    transform=axs[1].transAxes, fontsize=20, verticalalignment='top')
        
        if verbose:
            print(f"RMSE: {rmse:.3f}")
            print(f"Pearson correlation: {corr_coef:.3f} (p = {p_value:.3e})")

        if only_ranges:
            return {
            "key_type": key_type,
            "ranges": rangeList,
            "foundRange": foundRange,
            "RMSE": rmse,
            "corr_coef": corr_coef,
            "p_score": p_value,
            "slope": coef[0],
            "slope_err": slope_err
        }
        else: 
            plt.show()

            return {
            "key_type": key_type,
            "ranges": rangeList,
            "foundRange": foundRange,
            "fig": fig,
            "RMSE": rmse,
            "corr_coef": corr_coef,
            "p_score": p_value,
            "slope": coef[0],
            "slope_err": slope_err
            }


    return {
        "key_type": key_type,
        "ranges": rangeList,
        "foundRange": foundRange,
        "fig": fig
    }


def reco_cosine(dataPath, matrFile, xRange, rangesPSTAR, histName="hHitmap;1", acceptance=0.8, only_ranges=False, compare_pstar=True):
    """
    Perform cosine similarity analysis on measured ROOT histograms to detect peak ranges.

    This function processes all ROOT histogram files in a directory, computes the cosine similarity
    between each histogram and a reference matrix (system matrix), and determines the detected peak
    range for each measurement. Optionally, it plots the measured histograms and compares detected
    ranges to expected PSTAR ranges.

    Parameters
    ----------
    dataPath : str or Path
        Directory containing ROOT histogram files (.root).
    matrFile : np.ndarray
        System matrix used for cosine similarity calculation (columns correspond to reference vectors).
    xRange : np.ndarray
        Array of range values corresponding to bins in the system matrix and histograms.
    rangesPSTAR : dict
        Expected peak ranges (from PSTAR data), keyed by run ID or energy.
    histName : str, optional
        Name of the histogram inside the ROOT files to load. Default is "hHitmap;1".
    acceptance : float, optional
        Not currently used in the function. Reserved for future compatibility. Default is 0.8.
    only_ranges : bool, optional
        If True, skip plotting. Default is False.
    compare_pstar : bool, optional
        If True, compute metrics and plot detected vs. expected PSTAR ranges. Default is True.

    Returns
    -------
    rangeList : dict
        Detected peak positions keyed by run ID or energy.
    metrics : dict, optional
        Only returned if compare_pstar=True. Contains:
            - 'RMSE' : float
            - 'Correlation' : float
            - 'p_value' : float
            - 'slope' : float
    """

    dataPath = Path(dataPath)

    # Dictionary to store detected peak positions
    rangeList = {}
    key_type = None  # Will be set to "runid" or "energy" based on file naming

    # Prepare plotting axes if needed
    fig, axs = (None, None)
    if not only_ranges:
        fig, axs = plt.subplots(2, 1, figsize=(20, 15),
                                 gridspec_kw={'height_ratios': [2, 1]})
        cmap = plt.get_cmap("tab20")

    # --- Load histograms from ROOT files and process each ---
    for i, filePath in enumerate(sorted(dataPath.glob("*.root"))):
        parts = filePath.stem.split('_')

        # Skip background files
        if any("bg" in p.lower() for p in parts):
            continue

        # Determine key type based on filename
        if parts[0].startswith("run"):
            if key_type is None:
                key_type = "runid"
            key = int(parts[0].lstrip("run0"))
        else:
            if key_type is None:
                key_type = "energy"
            key = str(parts[0].rstrip("_"))

        if not filePath.is_file():
            continue

        try:
            # --- Load ROOT histogram ---
            simHist = raux.get_histo(path=filePath, histo_names=histName)
            image = simHist.vals.T[::-1].flatten()
    
            image /= image.sum()  # normalize histogram

     

            cosine_sim = (matrFile.T @ image) / (norm(image) * norm(matrFile, axis=0))
            cosine_sim_r = cosine_sim # reverse for plotting in correct range order

            # --- Detect peak range ---
            # --- Detect peak range and compatibility interval ---
            peak_index = np.argmax(cosine_sim_r)
            peak_value = cosine_sim_r[peak_index]

            rangeList[key] = xRange[peak_index]


        except Exception as e:
            print(f"Failed to load {filePath.name}: {e}")
            continue

        # --- Plot histogram cosine similarity if plotting enabled ---
        if not only_ranges:
            axs[0].plot(xRange, cosine_sim_r, label=f'{key} MeV', alpha=0.8, color=cmap(i % 20))
            axs[0].plot(xRange[peak_index], cosine_sim_r[peak_index], 'o', color='red')
            axs[0].set_xlim(xRange[0]-5, xRange[-1]+5)
            #axs[0].set_title("Measured Histograms (Cosine similarity)")
            axs[0].set_xlabel("Range [mm]")
            axs[0].set_ylabel("Cosine similarity")
            axs[0].legend()

    # --- Compare detected peaks with PSTAR expected values ---
    metrics = {}
    if compare_pstar:
        expected_ranges = []
        observed_ranges = []
        for run_id in rangeList:
            if run_id in rangesPSTAR:
                print(f'PSTAR = {rangesPSTAR[run_id]}, Peak = {rangeList[run_id]}')
                expected_ranges.append(rangesPSTAR[run_id])
                observed_ranges.append(rangeList[run_id])

        # Scatter plot of detected vs. expected
        axs[1].scatter(expected_ranges, observed_ranges, label='Data')

        # Linear fit for comparison
        coef, cov = np.polyfit(expected_ranges, observed_ranges, 1, cov=True)
        poly_fn = np.poly1d(coef)
        slope_err = np.sqrt(cov[0, 0])
        axs[1].plot(expected_ranges, poly_fn(expected_ranges),
                    label=f'slope = {coef[0]:.3f} ± {slope_err:.3f}', color='red')
        axs[1].legend()
        axs[1].set_ylabel('Observed peak [mm]')
        axs[1].set_xlabel('PSTAR expected [mm]')

        # Compute metrics
        expected = np.array(expected_ranges)
        observed = np.array(observed_ranges)
        y_pred = poly_fn(expected)
        rmse = np.sqrt(np.mean((observed - y_pred) ** 2))
        corr_coef, p_value = pearsonr(expected, observed)
        metrics['RMSE'] = rmse
        metrics['Correlation'] = corr_coef
        metrics['p_value'] = p_value
        metrics['slope'] = coef[0]

        # Display metrics on plot
        axs[1].text(0.05, 0.8, f"RMSE = {rmse:.4f}\nCorr. coeff = {corr_coef:.4f}",
                    transform=axs[1].transAxes, fontsize=20, verticalalignment='top')

        print(f"RMSE: {rmse:.3f}")
        print(f"Pearson correlation: {corr_coef:.3f} (p = {p_value:.3e})")

    if not only_ranges:
        plt.tight_layout()
        plt.show()

    # --- Return results ---
    if compare_pstar:
        return rangeList, metrics
    else:
        return rangeList


def create_system_matrix(dataPath, target_energies, histName="hHitmap;1", plot=True, interp=True, norm=True):
    """
    Create an interpolated and optionally normalized system matrix from ROOT hitmaps.

    This function loads ROOT histograms (hitmaps) from a folder, flattens them, optionally 
    normalizes the columns, interpolates the matrix to target energies, and optionally plots 
    the result.

    Parameters
    ----------
    dataPath : str or Path
        Path to folder containing ROOT hitmap files.
    target_energies : np.ndarray
        Energies at which to interpolate the system matrix.
    histName : str, optional
        Histogram name to extract from ROOT files. Default is "hHitmap;1".
    plot : bool, optional
        If True, display the system matrix as an image. Default is True.
    interp : bool, optional
        If True, interpolate the hitmaps to target_energies. Default is True.
    norm : bool, optional
        If True, normalize each column to sum to 1. Default is True.

    Returns
    -------
    sm_interp : np.ndarray
        Normalized and interpolated system matrix of shape (n_pixels, len(target_energies))
        if interp=True. Otherwise returns raw or normalized hitmaps.
    target_energies : np.ndarray
        Array of target energies corresponding to the columns of sm_interp. Only returned if interp=True.
    """
    all_hitmaps = []
    energies = []
    dataPath = Path(dataPath)
    for filePath in dataPath.glob("*.root"):
        if not filePath.is_file():
            continue

        # Read histogram
        simHist = raux.get_histo(path=filePath, histo_names=histName)
        hitmap = simHist.vals.T[::-1]  # shape (n_pixels_per_layer, n_layers)
        all_hitmaps.append(hitmap.flatten())

        # Extract energy from filename like 'hitmap_66dot76.root'
        energy_str = filePath.stem.split("_")[0]
        energy_value = float(energy_str.replace("dot", "."))
        energies.append(energy_value)

    all_hitmaps = np.array(all_hitmaps).T  # shape (n_pixels, n_files)
    energies = np.array(energies)

    # Interpolation
    interp_func = interp1d(energies, all_hitmaps, kind='cubic', axis=1, fill_value="extrapolate")
    hitmap_interpolated = interp_func(target_energies)

    # Normalize each column
    if norm:
        sm_interp = hitmap_interpolated / hitmap_interpolated.sum(axis=0)
    else:
        sm_interp = hitmap_interpolated

    # Optionally plot
    if plot:
        plt.figure(figsize=(12, 10))
        plt.imshow(sm_interp.T[::-1].T, aspect='auto', origin='lower',
                   extent=[target_energies[0], target_energies[-1], 0, sm_interp.shape[0]])
        plt.colorbar(label="Normalized intensity")
        plt.xlabel("Range [mm]")
        plt.ylabel("Detector pixel")
        plt.title("Interpolated System Matrix")

        # Add custom x-axis ticks
        num_ticks = 10
        ticks = np.linspace(target_energies[0], target_energies[-1], num_ticks)
        plt.xticks(ticks, [f"{t:.1f}" for t in ticks], rotation=45)
        plt.show()

    if interp:
        return sm_interp, target_energies
    
    elif not interp and norm:
        return all_hitmaps/all_hitmaps.sum(axis=0)
    else:
        return all_hitmaps

def score(result):
    """
    Compute a combined score for a single reconstruction result, used to rank iterations
    in a parameter scan.

    The score combines three metrics:
    - RMSE: lower is better
    - Correlation coefficient: closer to 1 is better
    - Slope: closer to 1 is better, weighted by slope error

    Parameters
    ----------
    result : dict
        Dictionary containing the reconstruction metrics. Expected keys:
        - 'RMSE' : float
        - 'corr_coef' : float
        - 'slope' : float
        - 'slope_err' : float

    Returns
    -------
    float
        Combined score where lower is better. Can be used to rank different iterations.
    """

    rmse = result["RMSE"]
    corr_diff = abs(1 - result["corr_coef"])
    slope_diff = abs(1 - result["slope"]) / (result["slope_err"] + 1e-6)

    return rmse + corr_diff + slope_diff


def find_best_params(dataPath, histName, matrFile, rangesPSTAR, xRange, nMin, nMax, nStep):
    """
    Perform a parameter scan over the number of iterations for MLEM reconstruction
    and identify the best iterations based on a scoring function.

    For each iteration count, calls `reco_mlem_from_sim` to compute metrics
    (RMSE, correlation coefficient, slope, etc.), computes a combined score, 
    and prints the top iterations ranked by score.

    Parameters
    ----------
    dataPath : str
        Path to the directory containing input ROOT hitmaps.
    histName : str
        Name of the histogram to process from the ROOT files.
    matrFile : str or np.ndarray
        System matrix used in MLEM reconstruction.
    rangesPSTAR : dict
        Optional dictionary of expected range values for comparison with reconstructed ranges.
    xRange : array-like
        Range of x-values for the reconstruction or fitting.
    nMin : int
        Minimum number of iterations to scan.
    nMax : int
        Maximum number of iterations to scan (exclusive).
    nStep : int
        Step size between iteration counts.

    Returns
    -------
    zip
        A zip object of top iteration numbers and their corresponding result dictionaries.

    Notes
    -----
    - Uses the `score()` function to rank iterations, combining RMSE, correlation coefficient, 
      and slope into a single scalar.
    - Top iterations are determined by sorting scores in ascending order (lower score = better).
    """
        
    results = {}
    for n in range(nMin, nMax, nStep):
        r = reco_mlem_from_sim(dataPath=dataPath,
                               histName=histName,
                               matrFile=matrFile,
                               rangesPSTAR=rangesPSTAR,
                               xRange=xRange,
                               nIter=n,
                               only_ranges=True)
        results[n] = r

    scores = {k: score(v) for k, v in results.items()}
    sorted_iters = sorted(scores.items(), key=lambda x: x[1])

    top_n = 10  
    best_iters = sorted_iters[:top_n]

    # extract iteration numbers and results
    best_iter_nums = [x[0] for x in best_iters]
    best_results = [results[i] for i in best_iter_nums]

    print("Top iterations:", best_iter_nums)
    for i, res in zip(best_iter_nums, best_results):
        print(f"Iteration {i}: RMSE={res['RMSE']:.5f}, Corr={res['corr_coef']:.5f}, Slope={res['slope']:.5f}")

    return(zip(best_iter_nums, best_results))