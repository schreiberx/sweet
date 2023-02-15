import os


def get_job_output_files(job_data):
    """
    Return a list of potential output files
    """

    if 'output.reference_filenames' in job_data:
        """
        We found some information in job_data itself and will use this
        """

        if job_data['output.reference_filenames'] != "":
            return job_data['output.reference_filenames'].split(";")

    """
    Try to find some output data files
    """
    if 'runtime.output_file_mode' in job_data:
        file_mode = job_data['runtime.output_file_mode']

        if file_mode in ["bin", ""]:
            ref_file_ending = ".sweet"
        elif file_mode == "csv":
            ref_file_ending = ".csv"
        else:
            raise Exception("Unknown output file mode")

    else:
        # Use binary one per default
        ref_file_ending = ".sweet"


    # Load all list of output files from reference job
    ref_files = []
    files = os.listdir(job_data['jobgeneration.job_dirpath'])
    for f in files:
        if f.endswith(ref_file_ending):
            ref_files.append(f)

    return ref_files

