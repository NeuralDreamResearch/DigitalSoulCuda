# DigitalSoul
## Unified platform for CPU, GPU, FPGA, Quantum Computing

  Digital Soul is a library, abstracting 4 distinict hardware platforms into one interface, thereby significantly reduce the effort need for developing comprehensive computational workflows. One key abstraction provides control over four different computational system, as of November 2023.


## Documentation

### Classes

#### 1. `LUTx_1` Class

##### Description

The `LUTx_1` class represents a LookUp Table (LUT) with a single output. It can be instantiated with a specified number of inputs and a logic ID. The class provides methods for computing the output on CPU or GPU, generating a unitary gate, and generating hardware description language (HDL) code.

##### Methods

- `computeCPU`: Computes the output using CPU.
- `computeGPU`: Computes the output using GPU (if available).
- `UnitaryGen`: Generates a unitary gate based on the LUT.
- `entityGen`: Generates HDL code for the LUT.

##### Example

```cpp
size_t inputs = 3;
size_t logicID = 5;
LUTx_1 lut(inputs, logicID);
bool inputValues[3] = {true, false, true};
bool output;
lut.computeCPU(inputValues, &output);
```

#### 2. `QN::Complex` Class

##### Description

The `QN::Complex` class represents a complex number. It provides methods for basic complex number operations such as addition, subtraction, multiplication, division, magnitude, argument, and conjugate.

##### Methods

- `operator+`, `operator-`, `operator*`, `operator/`: Basic complex number operations.
- `magnitude`: Computes the magnitude of the complex number.
- `arg`: Computes the argument of the complex number.
- `conj`: Computes the conjugate of the complex number.

##### Example

```cpp
QN::Complex<float> z1(1.0, 2.0);
QN::Complex<float> z2(2.0, -1.0);
QN::Complex<float> result = z1 + z2;
```

#### 3. `QN::Qudit` Class

##### Description

The `QN::Qudit` class represents a quantum qudit (quantum digit). It can be instantiated with a specified number of levels. The class provides methods for manipulating the state of the qudit.

##### Methods

- `oneHot`: Sets a specific level in the qudit state to 1.
- `loadStatevector`: Loads a new state vector into the qudit.
- `freeStatevector`: Frees the memory occupied by the state vector.
- `Psi`: Returns a pointer to the state vector.
- `numStates`: Returns the number of levels in the qudit.

##### Example

```cpp
size_t levels = 4;
Qudit<float> qudit(levels);
qudit.oneHot(2);
```

#### 4. `QN::Gate` Class

##### Description

The `QN::Gate` class represents a quantum gate. It can be instantiated with a specified dimension and an operator matrix. The class provides methods for loading a new operator matrix and transforming a quantum state using the gate.

##### Methods

- `loadOperator`: Loads a new operator matrix into the gate.
- `transform`: Transforms a given quantum state using the gate.

##### Example

```cpp
size_t dimension = 2;
Complex<float> operatorMatrix[4] = {Complex<float>(1, 0), Complex<float>(0, 0), Complex<float>(0, 0), Complex<float>(1, 0)};
Gate<float> quantumGate(dimension, operatorMatrix);
```

### Namespaces

#### 1. `QN` Namespace

The `QN` namespace encapsulates classes and functionality related to quantum computing.

##### Subnamespaces

- `metadata`: Contains predefined complex arrays representing Pauli matrices (X, Y, Z).


