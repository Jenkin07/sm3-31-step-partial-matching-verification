# SM3 Carry Absorption Probability Verification

Verification script for the carry absorption probability in the partial matching phase of the 31-step reduced-round SM3.

## Environment

- Windows 10/11, 64-bit
- Visual Studio 2022 (Desktop development with C++ workload)

## Build and Run

1. Open `SM3_Pr.vcxproj` in Visual Studio 2022.
2. Select **x64** as the target platform and **Release** (or **Debug**) as the configuration.
3. Build the project (Ctrl+Shift+B) and run it (Ctrl+F5).

## Notes

- The project is developed and tested with the Visual Studio 2022 IDE; no command-line build is required.
- The output prints the observed carry absorption probability for verification against the paper.
