# predictv30r
Predict models in R up to 3.2

This provides a functional interface to Predict R models so they may be served by OpenCPU.

## Public opencpu demo server
Uses the github webhook mechanism for updates as described here: https://www.opencpu.org/cloud.html
The name of the package in DESCRIPTION should match perfectly the name of the repository.
The github user triggering the hook should have a publicly available email address.

> For now, it seems the public cloud of opencpu doesn't accept the update of the apps via github webhook.
It asks for github CI - not clear whether that's webhook as it was 2 years ago and is still advertised online,
or whether it's github Action now.

The PREDICT3 package available on the opencpu cloud implements version 3.1 of the algorithm and can be used to test this version.
For now, version 3.2 tests can only be driven locally.

## Local opencpu instance
To drive the tests locally, clone the WintonCentre/predict-opencpu-server repository.
Open in Rstudio and run the lines from startup.R as appropriate.
If the URL is set to local, the tests should start running as soon as the REPL is connected to the test instance.

> Make sure to reload the environment fully and reimport the libraries one by one when you're hit with errors like "'' does not exist in the current directory".
