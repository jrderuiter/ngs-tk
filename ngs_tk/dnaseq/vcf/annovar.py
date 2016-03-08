

def get_annovar_fields(vcf_reader):
    """ Gets annovar fields in VCF from vcf_reader. """
    return [name for name, info in vcf_reader.infos.items()
            if 'ANNOVAR' in info.desc]


def summarize_annovar(rec, fields):
    """Summarizes annovar fields for record."""
    info = rec.INFO
    return {f: _flatten_value(info[f]) for f in fields}


def _flatten_value(val, delimiter=';'):
    if isinstance(val, (list, tuple)):
        return delimiter.join(filter(bool, val))
    else:
        return val
