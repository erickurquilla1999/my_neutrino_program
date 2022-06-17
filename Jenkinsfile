pipeline {
    triggers { pollSCM('') }  // Run tests whenever a new commit is detected.
    agent { dockerfile true } // Use the Dockerfile defined in the root Flash-X directory
    stages {

        //=============================//
    	// Set up submodules and amrex //
        //=============================//
    	stage('Prerequisites'){ steps{
		sh 'echo Hello World'
	}}

	//=========================//
	// Vacuum oscillation test //
	//=========================//
	stage('Vacuum Oscillation'){ steps{
	    dir('source'){
	        sh 'mkdir output'
		sh 'julia main.jl'
	    }
	    dir('tests'){
	        sh 'julia vacuum_oscillation_theorical_data.jl'
		sh 'python3 plot_vacuum_oscillation_test.py'
		archiveArtifacts artifacts: 'vacuum_test.pdf'
	    }
	}}

    } // stages{
} // pipeline{
