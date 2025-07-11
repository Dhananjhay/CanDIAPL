import csv

def read_csv(filepath):
    with open(filepath, newline='') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
    return rows


def generate_html_table(rows):
    header = rows[0]
    body_rows = rows[1:]

    header_html = "<tr>" + "".join(f"<th>{h}</th>" for h in header) + "</tr>"

    body_html = ""
    for row in body_rows:
        body_html += "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>"

    # Full table
    table_html = f"<table border='1' style='border-collapse: collapse;'>\n{header_html}\n{body_html}\n</table>"
    return table_html



def generate_html(table_path, output_html="viewer.html"):
    rows = read_csv(table_path)
    table_html = generate_html_table(rows)

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <title>Table Viewer</title>
        <style>
            th, td {{
                padding: 8px 12px;
                text-align: left;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
            }}
            thead tr {{
                position: sticky;
                top: 0;
                background: #f2f2f2;
            }}
        </style>
    </head>
    <body>
        {table_html}
    </body>
    </html>
    """

    with open(output_html, "w") as f:
        f.write(html_content)


table = snakemake.input["table"]
out_html = snakemake.output["viewer_html"]

generate_html(table_path=table, output_html=out_html)