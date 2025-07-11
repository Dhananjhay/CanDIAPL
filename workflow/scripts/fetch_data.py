import os
import pyvo

# TAP service endpoint
service = pyvo.dal.TAPService("https://gea.esac.esa.int/tap-server/tap")

# Load params
output_file = snakemake.output["out_path"]
ra = snakemake.params["ra"]
dec = snakemake.params["dec"]
radius = snakemake.params["radius"]

query = f"""
SELECT TOP 100
  source_id, ra, dec, parallax, phot_g_mean_mag,
  DISTANCE({ra}, {dec}, ra, dec) AS ang_sep
FROM gaiadr3.gaia_source
WHERE DISTANCE({ra}, {dec}, ra, dec) < {radius}
ORDER BY ang_sep ASC
"""

results = service.search(query)
table = results.to_table()

# Ensure output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Write to CSV
table.write(output_file, format="ascii.csv", overwrite=True)
