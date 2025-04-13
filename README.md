# predictv30r
Predict models in R up to 3.0

This provides a functional interface to Predict R models so they may be served by OpenCPU.

Uses the github webhook mechanism for updates as described here: https://www.opencpu.org/cloud.html
The name of the package in DESCRIPTION should match perfectly the name of the repository.
The github user triggering the hook should have a publicly available email address.

For now, it seems the public cloud of opencpu doesn't accept the update of the apps via github webhook.
It asks for github CI - not clear whether that's webhook as it was 2 years ago and is still advertised online,
or whether it's github Action now.

