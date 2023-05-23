import numpy as np
import matplotlib.pyplot as plt
import os

################################
########## filenames ###########
################################
def getFormattedDate(t, timescale):
    return "{:020.8f}".format(t * timescale);

def getFormattedIteration(niter):
    return "iter" + str(niter).zfill(3);

def getValFromVar(var, job, content):
   val = content;
   for i in range(len(var)):
       val = val[var[i]];
   return val[job];



################################
####### jobs organization ######
################################
def createGroups(list_jobs, groups_vars, jobs_data):

    groups = {};
    read_jobs = [];

    if groups_vars is None:
        groups[list_jobs[0]] = list_jobs
        return groups

    for job in list_jobs:

        new_group = True;

        if not (groups_vars is None):
            job_vars = [jobs_data[job][var] for var in groups_vars];
            for job2 in groups.keys():
                job2_vars = [jobs_data[job2][var] for var in groups_vars];
                if job_vars == job2_vars:
                    new_group = False;
                    break;

        if new_group:
            groups[job] = [job];
        else:
            groups[job2].append(job);

    return groups;

def getJobPath(data, job):
    return os.path.basename(data[job]["jobgeneration.p_job_dirpath"]);

###############
#### Plots ####
###############

## pretty names for labels
pretty_names = {
                   'runtime.timestep_size' : r'$\Delta t$',

                   'runtime.xbraid_max_levels' : r'$N_{\mathrm{levels}}$',
                   'runtime.xbraid_cfactor' : r'$m_c$',
                   'runtime.xbraid_nrelax' : r'$N_{\mathrm{relax}}$',
                   'runtime.xbraid_spatial_coarsening' : r'$M_{\mathrm{coarse}}$',
                   'runtime.xbraid_pt' : r'$N_{\mathrm{proc}}$',
                   'runtime.xbraid_viscosity_coefficient_fine' : r'$\nu_0$',
                   'runtime.xbraid_viscosity_coefficient_coarse' : r'$\nu_1$',
                   'runtime.xbraid_viscosity_coefficient_coarse1' : r'$\nu_1$',
                   'runtime.xbraid_viscosity_coefficient_coarse2' : r'$\nu_2$',
                   'runtime.xbraid_viscosity_coefficient_allcoarse' : r'$\nu_{\mathrm{coarse}}$',
                   'runtime.xbraid_viscosity_order_fine' : r'$q_0$',
                   'runtime.xbraid_viscosity_order_coarse' : r'$q_1$',
                   'runtime.xbraid_viscosity_order_coarse1' : r'$q_1$',
                   'runtime.xbraid_viscosity_order_coarse2' : r'$q_2$',

                   'err_L1' : r'$L_1$',
                   'err_L2' : r'$L_2$',
                   'err_L2_prog_phi' : r'$E_{\Phi, L_2}$',
                   'err_Linf' : r'$L_{\infty}$',
                   'rnorm_prog_phi_pert' : r'$E_{\Phi pert, R_{\mathrm{norm}}}$',
                   'rnorm_16_prog_phi_pert' : r'$E_{\Phi pert, R_{\mathrm{norm}} = 16}$',
                   'rnorm_32_prog_phi_pert' : r'$E_{\Phi pert, R_{\mathrm{norm}} = 32}$',
                   'rnorm_64_prog_phi_pert' : r'$E_{\Phi pert, R_{\mathrm{norm}} = 64}$',
                   'rnorm_128_prog_phi_pert' : r'$E_{\Phi pert, R_{\mathrm{norm}} = 128}$',
                   'rnorm_256_prog_phi_pert' : r'$E_{\Phi pert, R_{\mathrm{norm}} = 256}$',
                   'rnorm_prog_phi' : r'$E_{\Phi, R_{\mathrm{norm}}}$',
                   'rnorm_16_prog_phi' : r'$E_{\Phi, R_{\mathrm{norm}} = 16}$',
                   'rnorm_32_prog_phi' : r'$E_{\Phi, R_{\mathrm{norm}} = 32}$',
                   'rnorm_64_prog_phi' : r'$E_{\Phi, R_{\mathrm{norm}} = 64}$',
                   'rnorm_128_prog_phi' : r'$E_{\Phi, R_{\mathrm{norm}} = 128}$',
                   'rnorm_256_prog_phi' : r'$E_{\Phi, R_{\mathrm{norm}} = 256}$',
                   'rnorm_16_prog_vrt' : r'$E_{\xi, R_{\mathrm{norm}} = 16}$',
                   'rnorm_32_prog_vrt' : r'$E_{\xi, R_{\mathrm{norm}} = 32}$',
                   'rnorm_64_prog_vrt' : r'$E_{\xi, R_{\mathrm{norm}} = 64}$',
                   'rnorm_128_prog_vrt' : r'$E_{\xi, R_{\mathrm{norm}} = 128}$',
                   'rnorm_256_prog_vrt' : r'$E_{\xi, R_{\mathrm{norm}} = 256}$',

                   0: r'$0$',
                   10000.0 : r'$10^4$',
                   100000.0 : r'$10^5$',
                   1000000.0 : r'$10^6$',
                   '0': r'$0$',
                   '0.0': r'$0$',
                   '10000.0' : r'$10^4$',
                   '100000.0' : r'$10^5$',
                   '1000000.0' : r'$10^6$'
               };

short_names =  {
                   'runtime.timestep_size' : "dt",

                   'runtime.xbraid_max_levels' : "nlvl",
                   'runtime.xbraid_cfactor' : "cfact",
                   'runtime.xbraid_nrelax' : "nrlx",
                   'runtime.xbraid_spatial_coarsening' : "coars",
                   'runtime.xbraid_pt' : "pt",
                   'runtime.xbraid_viscosity_coefficient_fine' : "nu0",
                   'runtime.xbraid_viscosity_coefficient_coarse' : "nu1",
                   'runtime.xbraid_viscosity_coefficient_coarse1' : "nu1",
                   'runtime.xbraid_viscosity_coefficient_coarse2' : "nu2",
                   'runtime.xbraid_viscosity_coefficient_allcoarse' : "nucoars",
                   'runtime.xbraid_viscosity_order_fine' : "q0",
                   'runtime.xbraid_viscosity_order_coarse' : "q1",
                   'runtime.xbraid_viscosity_order_coarse2' : "q2",
               };

plot_attributes = {};
plot_attributes["color"] = plt.rcParams['axes.prop_cycle'].by_key()['color'];
plot_attributes["marker"] = ['+', 'o', 'x', 'D', 's', '^', 'v', '1', '2', '3', '4'];
plot_attributes["marker_thin"] = ['+', 'x', '1', '2', '3', '4'];
plot_attributes["linestyle"] = ['-', '--', '-.', ':'];

def getPlotAttribute(attribute, idx):
    return plot_attributes[attribute][idx % len(plot_attributes[attribute])];

def getPrettyName(v):

    ## special case: 1e+x values
    if type(v) is float:
        if v > 0:
            if np.log10(v) == int(np.log10(v)):
                return r'$10^{{{}}}$'.format(int(np.log10(v)))

    if v in pretty_names.keys():
        return pretty_names[v];
    else:
        return v;

def getShortName(v):

    if v in short_names.keys():
        return short_names[v];
    else:
        return v;


def getLabel(label_vars, job, jobs_data, xbraid_only_coarse = False):

    if label_vars is None:
        return ""

    if len(label_vars) > 1:
        label = "(";
        for var in label_vars:
            if "xbraid" in var and "," in str(jobs_data[job][var]) and xbraid_only_coarse:
                label += str(getPrettyName(jobs_data[job][var].split(",")[-1])) + ",";
            else:
                label += str(getPrettyName(jobs_data[job][var])) + ",";
        label = label[:-1] + ")";
    else:
        var = label_vars[0];
        label = getPrettyName(var) + r'$ = $' + str(getPrettyName(jobs_data[job][var]));
    return label;

def setFilenameFilterVars(figname, filter_vars, filter_vals):

    figname = figname + "_";

    ivar = 0;
    for var in filter_vars:

        if hasattr(filter_vals[ivar], "__len__") and not isinstance(filter_vals[ivar], str):
            if len(filter_vars) < 5:
                figname += var
            else:
                figname += getShortName(var)
            for val in filter_vals[ivar]:
                figname += "_" + str(val);
            figname += "__"
        else:
            if len(filter_vars) < 5:
                figname += var;
            else:
                figname += getShortName(var)
            figname += "_" + str(filter_vals[ivar]) + "__";
        ivar += 1;

    return figname;


def setTitleFilterVars(title, filter_vars, filter_vals):

    title += "\n";

    ivar = 0;
    for var in filter_vars:

        if ivar > 0 and ivar % 6 == 0:
            title += "\n";

        title += getPrettyName(var) + "=" + getPrettyName(str(filter_vals[ivar]), var_principal = var[-1]) + "; ";

        ivar += 1;

    var = var[:-2];

    return title;

def setTitleLegendVars(title, legend_vars):

    title += "\n";
    title += "In function of (";

    for var in legend_vars:
        title += getPrettyName(var) + ",";
    title = title[:-1] + ")";

    return title;

def setTitleLegend(legend_vars):

    if len(legend_vars) == 1:
         return ""

    title = "(";

    for var in legend_vars:
        title += getPrettyName(var) + ",";
    title = title[:-1] + ")";

    return title;



def filterJob(job, filter_vars, filter_vals, jobs_data, keep):

    ok = True;

    ivar = 0;
    for var in filter_vars:
        if hasattr(filter_vals[ivar], "__len__") and not isinstance(filter_vals[ivar], str):
            if not (jobs_data[job][var] in filter_vals[ivar]):
                ok = False;
                break;
        else:
            if not (jobs_data[job][var] == filter_vals[ivar]):
                ok = False;
                break;

        ivar += 1;

    if keep:
        return ok;
    else:
        return not ok;

def filterJobs(list_jobs, filter_vars, filter_vals, jobs_data, keep = True):

    if filter_vars is None:
        return list_jobs

    list_jobs_copy = list_jobs.copy();
    for job in list_jobs_copy:
        if not filterJob(job, filter_vars, filter_vals, jobs_data, keep):
            list_jobs.remove(job);

    return list_jobs;


def excludeJobs(list_jobs, exclude_vars, all_exclude_vals, jobs_data):
    for exclude_vals in all_exclude_vals:
        list_jobs = filterJobs(list_jobs, exclude_vars, exclude_vals, jobs_data, keep = False)
    return list_jobs

def sortJobs(list_jobs, sort_vars, jobs_data):

    if sort_vars is None:
        return list_jobs

    sorted_jobs = [];
    for job in list_jobs:
        arr = [];
        for var in sort_vars:
            arr.append(jobs_data[job][var]);
        arr.append(job);
        sorted_jobs.append(arr);
    sorted_jobs.sort(key = lambda x: x[:]);
    sorted_jobs = [row[len(sort_vars)] for row in sorted_jobs];

    return sorted_jobs



###################
## SPECTRAL DATA ##
###################

## get index of mode (n, m)
def getArrayIndexByModes(n, m, N_max):

    assert n >= 0;
    assert n >= m;

    return (m * (2 * N_max - m + 1) >> 1)  + n;

## get max abs considering the first rnorm modes
def getMaxAbsRnorm(u, rnorm, N_max, verbose = False):

    err = 0;

    ##########for m in range(int(rnorm)):
    ##########    for n in range(m, int(rnorm)):
    ##########        idx = getArrayIndexByModes(n, m, N_max);
    ##########        err = np.max([err, np.abs(u[idx] * np.conj(u[idx]))]);
    ##########        if verbose:
    ##########            print(m, n, idx, u[idx], err);

    if u is None:
        return None

    for m in range(int(rnorm)):
        for n in range(m, int(rnorm)):
            if m < u.shape[0] and n < u.shape[1]:
                err = np.max([err, np.abs(u[m, n] * np.conj(u[m, n]))]);

    return float(err);

##idem but also return index
def getMaxAbsRnormWithIndex(u, rnorm, N_max):

    err = 0;
    i = -1;

    for m in range(int(rnorm)):
        for n in range(m, int(rnorm)):
            idx = getArrayIndexByModes(n, m, N_max);
            err = np.max([err, np.abs(u[idx] * np.conj(u[idx]))]);
            if np.abs(err - np.abs(u[idx] * np.conj(u[idx]))) < 1e-13:
                i = idx;

    return float(err), i;

