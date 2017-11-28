README for nais_processer.py
----------------------------

1. See that you have all the python-packages installed.

    sudo apt-get install python-matplotlib python-scipy python-numpy python-pandas python-easygui

2. Navigate to the program directory and run the program

    cd /path/to/nais_processer/
    python nais_processer.py

3. The program asks for the following information:

    Path to files
    Start date (yyyymmdd)
    End date (yyyymmdd)
    Inverter resolution (lores/hires)
    Sampleflow (lpm)
    STP correction (yes/no)
    Tubeloss correction (yes/no)
    Inlet length (m)
    Inlet diameter (m)
    Delimiter in files (\t/,/etc.)
    Date format in filenames (e.g. %Y%m%d or %Y-%m-%d)
    Records/diagnostics filename ending
    Ion filename ending
    Particle filename ending
    Temperature/pressure recorded? (yes/no)
    Deviation from UTC in hours (negative or positive integer)

    NOTE:
      - Models newer than NAIS12 record temperature/pressure and have a 54 lpm sampleflow. 
      - Older models and AIS have 60 lpm sampleflow and do not record temperature/pressure.
      - For older models you can define your own temperature/pressure in the code (default: 1 atm, 25 C)
      - Find out the inverter resolution by for example counting the mobility bins in the ion data file (lores has 28 bins and hires has 109 bins).
      - STP correction = convert the data to standard conditions (1 atm, 25 C)
      - Tubeloss correction = correct for laminar diffusion losses to walls of the inlet tube
      - Check carefully for typos, just basic python error handling.

    When done press OK.

4. Choose time from the list of parameters (use begin time) and press OK

5. If temperature/pressure recorded? = 'yes' you will also be asked for temperature and pressure respectively.

6. The program creates two folders in the "Path to files"-folder called: sum and figures. 

    sum-folder contains the number-size distributions as sum-matrices.

    sum-matrix explanations (python indexing):

        sum_matrix[1:,0] = matlab datenum
        sum_matrix[0,2:] = bin geometric mean diameters [m]
        sum_matrix[1:,2:] = normalized concentrations, dNdlogDp [cm-3]
        sum_matrix[1:,2] = Integrated total number concentration [cm-3]

    figures-folder contains the surface plots of the number-size distributions.
