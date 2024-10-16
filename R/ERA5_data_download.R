test_cds_api <- function() {

  ###
  # TODO: This doesn't work yet, in R, but does work in terminal, should fix later!
  ###

  python_code <- '
import cdsapi

def test_connection():
    c = cdsapi.Client()
    try:
        c.retrieve(
            "reanalysis-era5-single-levels",
            {
                "product_type": "reanalysis",
                "format": "netcdf",
                "variable": "2m_temperature",
                "year": "2022",
                "month": "01",
                "day": "01",
                "time": "12:00",
                "area": [90, -180, -90, 180],
            },
            "test_download.nc"
        )
        print("CDS API connection successful!")
        return True
    except Exception as e:
        print(f"Error: {str(e)}")
        return False
'

  reticulate::py_run_string(python_code)
  reticulate::py$test_connection()
}

# Create a function to run the Python code
download_era5_data <- function(start_year, end_year) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    install.packages("reticulate")
  }

  # First, test the API connection
  if (test_cds_api()) {
    # If successful, proceed with the data download
    # Ensure the cdsapi Python package is installed
    # reticulate::py_install("cdsapi")

    # Python code as a string
    python_code <- '
import cdsapi

def download_data(start_year, end_year):
    c = cdsapi.Client()

    for year in range(start_year, end_year + 1):
        c.retrieve(
            "reanalysis-era5-single-levels",
            {
                "product_type": "reanalysis",
                "format": "netcdf",
                "variable": ["2m_temperature", "relative_humidity", "leaf_area_index_high_vegetation", "surface_solar_radiation_downwards", "volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2"],
                "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"],
                "year": year,
                "day": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31"],
                "time": ["00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"],
                "area": [61.9, 24.2, 61.8, 25.3],
            },
            f"{year}download.nc"
        )

        c.retrieve(
            "reanalysis-era5-single-levels",
            {
                "product_type": "reanalysis",
                "format": "netcdf",
                "variable": "2m_dewpoint_temperature",
                "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"],
                "year": year,
                "day": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31"],
                "time": ["00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"],
                "area": [61.9, 24.2, 61.8, 25.3],
            },
            f"{year}download_td.nc"
        )
'

    # Use reticulate to run the Python code
    reticulate::py_run_string(python_code)

    # Call the Python function with R arguments
    reticulate::py$download_data(as.integer(start_year), as.integer(end_year))
  } else {
    print("CDS API connection failed. Please check your API key and internet connection.")
  }
}
