
pipeline {
    agent any

    stages {
        stage('Test') {
            steps {
                echo 'Testing..'
                checkout scm
		sh 'echo "hello world"'
		sh 'source ./activate.sh'
		sh 'export'
            }
        }
    }
}


