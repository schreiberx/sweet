import os
import sys
import glob





def add_source_files(env, p):
    #
    # Add an action to move any module files
    #
    def moveModFiles(target=None, source=None, env=None):
        import glob, os, os.path
        targetdir = target[0].dir
        for t in target:
            if t.name[-4:] == '.mod':
                os.rename(t.name,os.path.join(str(targetdir),t.name))

    co = p.get_program_specific_options()

    if co == None:
        return

    sweet_root = env['MULE_SOFTWARE_ROOT']+'/'

    #
    # Load and iterate over all source code files and directories
    #
    fad_list = co['compile_files_and_dirs']
    for k in fad_list:
        abs_k = sweet_root+'/'+k
        if os.path.isdir(abs_k):
            print("Processing additional directory '"+abs_k+"'")

            cpp_files = glob.glob(abs_k+'/*.cpp')
            for i in cpp_files:
                print(" + Adding source file "+i)
                filerelpath = i.replace(sweet_root+"/", '')

                # SWE REXI special file handling for threaded parallelization over the REXI sum
                filetmp = os.path.basename(filerelpath)
                if 'l_rexi' in filetmp or 'lg_rexi' in filetmp or 'lc_rexi' in filetmp:
                    if p.rexi_thread_parallel_sum=='enable':
                        env_omp = env.Clone()
                        env_omp.Append(CXXFLAGS = ' -fopenmp')
                        env_omp.src_files.append(env_omp.Object(filerelpath))
                    else:
                        print(filerelpath)
                        env.src_files.append(env.Object(filerelpath))
                else:
                    env.src_files.append(env.Object(filerelpath))

            fortran_files = glob.glob(abs_k+'/*.f90')
            for i in fortran_files:
                filerelpath = i.replace(sweet_root+'/', '')

                obj = env.Object(filerelpath)
                env.src_files.append(obj)
                #env.AddPostAction(obj, moveModFiles)

        elif os.path.isfile(abs_k):
            print("Processing additional file '"+abs_k+"'")

            i = abs_k
            print(" + Adding source file "+i)
            filerelpath = i.replace(sweet_root+"/", '')

            obj = env.Object(filerelpath)
            env.src_files.append(obj)

        else:
            raise Exception("Error file processing file or directory '"+abs_k+"'")


Import('env', 'p')
add_source_files(env, p)
Export('env', 'p')

