
pipeline {
    agent any

    stages {
	parallel {
		stage('Test 1') {
		    steps {
			checkout scm
			sh '/bin/bash -c "source ./activate.sh; ./tests/05_jobs_run_directly_compile_error/test.py"'
		    }
		}

		stage('Test 2') {
		    steps {
			checkout scm
			sh '/bin/bash -c "source ./activate.sh; ./tests/05_jobs_run_directly_compile_error/test.py"'
		    }
		}
	}
    }
}


