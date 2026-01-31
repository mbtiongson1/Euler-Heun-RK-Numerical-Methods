# utils.py

def print_table(headers, rows):
    """
    Print a table nicely formatted for console output.
    Floats are truncated to 6 decimal places.
    """
    print(" | ".join(headers))
    print("-" * (len(headers) * 12))
    for row in rows:
        formatted_row = []
        for val in row:
            if isinstance(val, float):
                formatted_row.append(f"{val:.6f}")
            else:
                formatted_row.append(str(val))
        print(" | ".join(formatted_row))


def print_table_csv(headers, rows):
    """
    Print a table in CSV format (comma-separated),
    ready for Excel copy-paste.
    Floats formatted to 6 decimals, percentages preserved as strings.
    """
    print(",".join(headers))
    for row in rows:
        formatted_row = []
        for val in row:
            if isinstance(val, str) and val.endswith('%'):
                formatted_row.append(val)
            elif isinstance(val, float):
                formatted_row.append(f"{val:.6f}")
            else:
                formatted_row.append(str(val))
        print(",".join(formatted_row))