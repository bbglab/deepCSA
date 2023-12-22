

def add_filter(old_filt, add_filt, filt_name):
    """
    old filt is the current FILTER field value
    add_filt is a boolean, either True or False
    filt_name is the name that should be added in the FILTER field in case the add_filt value is True
    """
    if add_filt:
        if old_filt == "PASS":
            return filt_name
        old_filt += ";" + filt_name

    return ";".join( sorted(old_filt.split(";")) )

def to_int_if_possible(string):
    try:
        int(string)
        return int(string)
    except ValueError:
        return None
