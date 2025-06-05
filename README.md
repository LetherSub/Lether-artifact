# Lether-artifact

This code is based on the [LaZer library](https://github.com/lazer-crypto/lazer), so the dependencies are the same as those required by LaZer.

## Build Instructions

Follow these steps to compile and run the Lether demo:

```bash
# Clean and compile the root project
make clean
make all -j8

# Enter the Lether demo directory
cd ./demos/Lether

# Clean and compile the Lether demo
make clean
make

# Run the demo
./Lether-demo