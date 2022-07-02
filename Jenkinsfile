
pipeline {
    agent any

    stages {
        stage('Test') {
            steps {
                checkout scm
		sh 'source ./activate.sh; ./tests/05_jobs_run_directly_compile_error/test.py'
            }
        }
    }
}


