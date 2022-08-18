import numpy as np
import struct
import os

def read_bin_file(filename):
    f = open(filename, 'rb');

    ## LINE 1
    c = np.fromfile(f, dtype = '|S6', count = 1).astype("|U6")
    assert c == "SWEET\n";

    ## LINE 2
    skip = 6;
    f.seek(skip, os.SEEK_SET);
    c = np.fromfile(f, dtype = '|S18', count = 1).astype("|U18")
    assert c == "DATA_TYPE SH_DATA\n";

    ## LINE 3
    skip += 18;
    f.seek(skip, os.SEEK_SET);
    c = np.fromfile(f, dtype = '|S12', count = 1).astype("|U12")
    assert c == "MODES_M_MAX ";
    modes_m_max = [];
    s = ''
    cnt = 0;
    while True:
        s = np.fromfile(f, dtype = '|S1', count = 1).astype("|U1")
        cnt += 1;
        if s == '\n':
            break;
        modes_m_max.append(s);
    m_max = '';
    for s in modes_m_max:
        m_max = m_max + str(s[0]);
    m_max = int(m_max)

    ## LINE 4
    skip += 12 + cnt;
    f.seek(skip, os.SEEK_SET);
    c = np.fromfile(f, dtype = '|S12', count = 1).astype("|U12")
    assert c == "MODES_N_MAX ";
    modes_n_max = [];
    s = ''
    cnt = 0;
    while True:
        s = np.fromfile(f, dtype = '|S1', count = 1).astype("|U1")
        cnt += 1;
        if s == '\n':
            break;
        modes_n_max.append(s);
    n_max = '';
    for s in modes_n_max:
        n_max = n_max + str(s[0]);
    n_max = int(n_max)

    ## LINE 5
    skip += 12 + cnt;
    f.seek(skip, os.SEEK_SET);
    c = np.fromfile(f, dtype = '|S10', count = 1).astype("|U10")
    assert c == "GRID_TYPE ";
    grid_type = [];
    s = ''
    cnt = 0;
    while True:
        s = np.fromfile(f, dtype = '|S1', count = 1).astype("|U1")
        cnt += 1;
        if s == '\n':
            break;
        grid_type.append(s);
    gt = '';
    for s in grid_type:
        gt = gt + str(s[0]);

    ## LINE 6
    skip += 10 + cnt;
    f.seek(skip, os.SEEK_SET);
    c = np.fromfile(f, dtype = '|S13', count = 1).astype("|U13")
    assert c == "NUM_ELEMENTS ";
    num_elements = [];
    s = ''
    cnt = 0;
    while True:
        s = np.fromfile(f, dtype = '|S1', count = 1).astype("|U1")
        cnt += 1;
        if s == '\n':
            break;
        num_elements.append(s);
    ne = '';
    for s in num_elements:
        ne = ne + str(s[0]);

    ## LINE 7
    skip += 13 + cnt;
    f.seek(skip, os.SEEK_SET);
    c = np.fromfile(f, dtype = '|S4', count = 1).astype("|U4")
    assert c == "FIN\n";

    ## LINE 8+
    skip += 4;
    f.seek(skip, os.SEEK_SET);
    data  = np.fromfile(f, dtype = np.complex_)

    return data, m_max, n_max;

###read_file(filename)

####SWEET
####DATA_TYPE SH_DATA
####MODES_M_MAX 31
####MODES_N_MAX 31
####GRID_TYPE GAUSSIAN
####NUM_ELEMENTS 528
####FIN
