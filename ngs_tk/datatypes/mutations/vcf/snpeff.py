
def get_snpeff_fields(reader):
    field_str = reader.infos['ANN'].desc.split('\'')[1]
    return field_str.split(' | ')

def summarize_snpeff(record, fields):
    # Check if we can have multiple entries.

    if 'ANN' not in record.INFO:
        return None

    entries = (dict(zip(fields, value_str.split('|')))
               for value_str in record.INFO['ANN'])

    return tuple(entries)
