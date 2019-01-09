## WCSim Analysis for ANNIE
In order to keep analysis code manageable this repository provides an analysis class for performing analysis of the outputs of WCSim. 

## Running the code
1. Get yourself a copy of the repository code. 

2. ROOT6 (which the code relies on) requires the WCSim headers be directly `#included` in the analysis code, so if you don't have a copy of ANNIE's WCSim repository, you should get one, or at the very least the `include` directory, the `libWCSimRoot.so` library and the `WCSimRootDict_rdict.pcm` pre-compiled module. 

3. Since these paths are hard-coded into the analsis code, the required directory structure is: 
    * wcsim     [source code directory - required]
        * include [headers - required]
        * src     [source files - optional]
    * build     [build directory - optional]
    * analysis   [analysis code - required]

4. Move to the analysis folder
```cd analysis```

5. Start ROOT6
```$ root -l```

6. Load the `analysiscaller.cxx` file
```root [0] .L /path/to/analysis/analysiscaller.cxx+g```

7. Instantiate a _WCSimAnalysis_
```root [1] WCSimAnalysis* theana = analysiscaller(/path/to/files)```
Where the path points to a directory containing wcsim_*.root files. All files in this directory and any subdirectories will be added to the analysis TChain.

8. This performs the analaysis as `DoAnalysis()` is called in the constructor. 

9. Inspect the output of the analysis.

10. Delete the analysis class when done
`root [2] delete theana`
[note this will generate a segfault if the canvases are closed before deleting theana].

## Main Files and Analysis Structure

#### Intialisation and Utility

* **`wcsimanalysis.hh`** defines the main class _WCSimAnalysis_, which contains the analysis functions and a collection of members useful for passing between the various functions. The function definitions are in the files `#included` at the end of the file.

* **`analysiscaller.cxx`** is the top level analysis file. It's main job is to instantiate a _WCSimAnalysis_ and call it's `DoAnalysis()` function, which is the main top level analysis function.

* **`doanalysis.cxx`** defines that function. It's job is to call in turn the necessary steps to peform the analysis.

* **`utilityfuncs.cxx`** contains definitions for the first few functions called. These are general file- and environment-level functions, such as loading the `libWCSimRoot.so` library, creating a `TChain` and adding the input files, loading TChain entries with proper handling file transitions, and so on.

* **`makepmtmaps.cxx`** defines a single function that produces a set of std::map's to convert a PMT ID into an (x,y) pair that represents the PMT's position. These maps are used, for example, to generate a 2D map of hit frequency on the tank wall and caps _a la_ 'SK style'. 

* **`[element]hists.cxx`** defines histograms, their fill and draw calls - split by element (tank/mrd/veto). Define all your histograms and canvases here, and add the pointers to the WCSimAnalysis class. 

#### Looping over Events
The main analysis consists of functions to call before, during, and after the event loop. These `Do[Element][Stage]` functions themselves contain a collection of _function calls_, with the analysis proper being contained within those called functions. Each analsis can thus be enabled or disabled just by commenting out the relevant function calls. 

The **`DoXXXPreLoop`** and **`DoXXXPostLoop`** functions are called before and after the event loop, respectively.

#### Looping over Hits and Digits
During the loop over the events five functions are called for each detector element. 

* **`DoXXXEventWide`** performs any analysis that just uses top level information - such as the event header, counts of hits, digits, hit tubes etc.

* **`DoXXXPreHitLoop`** is next called before looping over the hits and digits. This is for doing initialisation steps necessary prior to the loop over hits/digits. 

* **`DoXXXTrueHits`** performs a loop over the true hits. It retrieves the hit and all associated information, and passes it to calls to the analysis functions that analyze hits.

* **`DoXXXDigitHits`** similarly performs a loop over digit hits; each time retrieving the digi and passing it to functions that analyze digits. 

* **`DoXXXPostHitLoop`** is called at the end of the loop over hits and digits. Any analysis that requires collections filled during the loops is performed here.

* In the future there may also be a `DoXXXTracks` function to loop over tracks, - so far no analysis is done with the tracks collection. 

#### LAPPDs

This structure has no present handles for LAPPD information; this information is stored in separate files with a different structure, so it's implementation requires consideration.

____

## Contributing to Development
Moving forward with the analysis means creating new function calls in `Do----` functions.
