# KyberUPKE
Updatable public key encryption from Kyber

The original code is from https://github.com/GiacomoPope/kyber-py. We fixed the rho for matrix A and then added a function unKey in kyber.py, which allowed updating the public key and private key according to a random 32-byte sigma.

This project is only used to test the UPKE scheme from Kyber and cannot be applied to the real implementations. 
