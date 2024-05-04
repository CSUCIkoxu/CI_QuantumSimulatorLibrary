using System;
using System.Numerics;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;    //Need to import Mathnet into csharp.proj using nuGet or any other method
using MathNet.Numerics.LinearAlgebra.Complex;

namespace qSim {
    public class qSimulator {
        public static Complex[][,] operatorMat = new Complex[11][,];
            // [[[new Complex(1/Math.Sqrt(2),0),new Complex(1/Math.Sqrt(2),0)],[new Complex(1/Math.Sqrt(2),0),new Complex(-1/Math.Sqrt(2),0)]],  //H
            // [[new Complex(1,0),new Complex(0,0)],[new Complex(0,0),new Complex(1,0)]],  //I
            // [[new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)],[new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0)],[new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0)]],  //CNOT
            // [[new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)],[new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0)],[new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0)]],  //CCNOT
            // [[new Complex(1,0),new Complex(0,0)],[new Complex(0,0),new Complex(1,0)]],   //R STUB
            // [[new Complex(1,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,1)]],  //S or P
            // [[new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0)],[new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(1,0)]],  //SWAP
            // [[new Complex(1,0),new Complex(0,0)],[new Complex(0,0),new Complex(1/Math.Sqrt(2),1/Math.Sqrt(2))]],    //T
            // [[new Complex(0,0),new Complex(1,0)],[new Complex(1,0),new Complex(0,0)]],  //X
            // [[new Complex(0,0),new Complex(0,-1)],[new Complex(0,1),new Complex(0,0)]], //Y
            // [[new Complex(1,0),new Complex(0,0)],[new Complex(0,0),new Complex(-1,0)]]]; //Z
        public Matrix<Complex>[] operators;

        private MathNet.Numerics.LinearAlgebra.Vector<Complex> stateVector;
        private int qbitQuantity;
        private Random rand;

        //Constructors//////////
        public qSimulator() 
            : this(4)
        {}
        public qSimulator(int qBitNum, int seed = -1) {
            this.stateVector = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense((int) Math.Pow(2,qBitNum));
            this.stateVector[0] = 1;
            this.qbitQuantity = qBitNum;

            operatorMat[0] = new Complex[,] {{new Complex(1/Math.Sqrt(2),0),new Complex(1/Math.Sqrt(2),0)},{new Complex(1/Math.Sqrt(2),0),new Complex(-1/Math.Sqrt(2),0)}};  //H
            operatorMat[1] = new Complex[,] {{new Complex(1,0),new Complex(0,0)},{new Complex(0,0),new Complex(1,0)}};  //I
            operatorMat[2] = new Complex[,] {{new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)},{new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0)},{new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0)}};  //CNOT
            operatorMat[3] = new Complex[,] {{new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)},{new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0)},{new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0)}};  //CCNOT
            operatorMat[4] = new Complex[,] {{new Complex(1,0),new Complex(0,0)},{new Complex(0,0),new Complex(1,0)}};   //R STUB
            operatorMat[5] = new Complex[,] {{new Complex(1,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,1)}};  //S or P
            operatorMat[6] = new Complex[,] {{new Complex(1,0),new Complex(0,0),new Complex(0,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(0,0)},{new Complex(0,0),new Complex(0,0),new Complex(1,0),new Complex(1,0)}};  //SWAP
            operatorMat[7] = new Complex[,] {{new Complex(1,0),new Complex(0,0)},{new Complex(0,0),new Complex(1/Math.Sqrt(2),1/Math.Sqrt(2))}};    //T
            operatorMat[8] = new Complex[,] {{new Complex(0,0),new Complex(1,0)},{new Complex(1,0),new Complex(0,0)}};  //X
            operatorMat[9] = new Complex[,] {{new Complex(0,0),new Complex(0,-1)},{new Complex(0,1),new Complex(0,0)}}; //Y
            operatorMat[10] = new Complex[,] {{new Complex(1,0),new Complex(0,0)},{new Complex(0,0),new Complex(-1,0)}}; //Z

            this.operators = new Matrix<Complex>[operatorMat.Length];
            for (int i = 0; i < operatorMat.Length; i++) {
                this.operators[i] = Matrix<Complex>.Build.DenseOfArray(operatorMat[i]);
            }

            if (seed != -1) {
                this.rand = new Random(seed);
            } else {
                this.rand = new Random();
            }
        }

        //Public Functions//////////
        
        //Warning: Deletes any operations you may have performed before hand
        public void newShape(int newqBitNum) {
            MathNet.Numerics.LinearAlgebra.Vector<Complex> newStateVector = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense((int) Math.Pow(2,newqBitNum));
            this.stateVector = newStateVector;
        }
        public void runOp(opCodes opNum, int[] bits, double[] parameters){
            Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});
            switch (opNum) {
                case opCodes.M:
                    this.stateVector = measureBit(this.stateVector, bits[0]);
                    break;
                case opCodes.H:
                    //Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});

                    for (int i = 0; i < this.qbitQuantity; i++) {
                        if (i == bits[0]) {
                            currOp = tensorProd(currOp, this.operators[0]); //H
                        } else {
                            currOp = tensorProd(currOp, this.operators[1]); //I
                        }
                    }
                    
                    currOp.Multiply(this.stateVector, this.stateVector);
                    break;
                case opCodes.I:
                    break;
                case opCodes.CNOT:
                    //Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});

                    for (int i = 0; i < this.qbitQuantity; i++) {
                        if (i == bits[0]) {
                            currOp = tensorProd(currOp, this.operators[2]); //CNOT
                        } else {
                            currOp = tensorProd(currOp, this.operators[1]); //I
                        }
                    }

                    MathNet.Numerics.LinearAlgebra.Vector<Complex> tempStateVec = this.stateVector.Clone();
                    int control = bits[0];
                    int target = bits[1];

                    tempStateVec = swapBits(tempStateVec, target, control + 1); //Check - possible logic error
                    tempStateVec = currOp.Multiply(tempStateVec);
                    tempStateVec = swapBits(tempStateVec, control + 1, target);

                    this.stateVector = tempStateVec;
                    break;
                case opCodes.CCNOT:
                    break;
                case opCodes.R:
                    //Determine the Axis of rotation
                    //parameter: [pauli axis, rotation(degrees)]
                    //Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});
                    Matrix<Complex> currAxis;
                    switch((PauliAxis) ((int) parameters[0])) {
                        case PauliAxis.X:
                            currAxis = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{new Complex(Math.Cos(parameters[1]/2),0), new Complex(0, -1 * Math.Sin(parameters[1]/2))},
                                                                        {new Complex(0, -1 * Math.Sin(parameters[1]/2)), new Complex(Math.Cos(parameters[1]/2),0)}});
                            break;
                        case PauliAxis.Y:
                            currAxis = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{new Complex(Math.Cos(parameters[1]/2),0), new Complex(-1 * Math.Sin(parameters[1]/2), 0)},
                                                                        {new Complex(Math.Sin(parameters[1]/2), 0), new Complex(Math.Cos(parameters[1]/2),0)}});
                            break;
                        case PauliAxis.Z:
                            currAxis = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{new Complex(Math.Cos(-1* parameters[1]/2),Math.Sin(-1* parameters[1]/2)), new Complex(0, 0)},
                                                                        {new Complex(0, 0), new Complex(Math.Cos(parameters[1]/2),Math.Cos(parameters[1]/2))}});
                            break;
                        default:
                            currAxis = this.operators[1]; //I
                            break;
                    }

                    //Build Operator Matrix
                    for (int i = 0; i < this.qbitQuantity; i++) {
                        if (i == bits[0]) {
                            currOp = tensorProd(currOp, currAxis); //Based on selected Pauli Axis and theta
                        } else {
                            currOp = tensorProd(currOp, this.operators[1]); //I
                        }
                    }
                    
                    //Apply Operator
                    currOp.Multiply(this.stateVector, this.stateVector);
                    break;
                case opCodes.S:
                    //Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});

                    for (int i = 0; i < this.qbitQuantity; i++) {
                        if (i == bits[0]) {
                            currOp = tensorProd(currOp, this.operators[5]); //S
                        } else {
                            currOp = tensorProd(currOp, this.operators[1]); //I
                        }
                    }
                    
                    currOp.Multiply(this.stateVector, this.stateVector);
                    break;
                case opCodes.SWAP:
                    //Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});

                    for (int i = 0; i < this.qbitQuantity; i++) {
                        if (i == bits[0]) {
                            currOp = tensorProd(currOp, this.operators[6]); //SWAP
                        } else {
                            currOp = tensorProd(currOp, this.operators[1]); //I
                        }
                    }

                    MathNet.Numerics.LinearAlgebra.Vector<Complex> tempStateVecS = this.stateVector.Clone();
                    int controlS = bits[0];
                    int targetS = bits[1];

                    tempStateVecS = swapBits(tempStateVecS, targetS, controlS + 1); //Check - possible logic error
                    tempStateVecS = currOp.Multiply(tempStateVecS);
                    tempStateVecS = swapBits(tempStateVecS, controlS + 1, targetS);

                    this.stateVector = tempStateVecS;
                    break;
                case opCodes.T:
                    //Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});

                    for (int i = 0; i < this.qbitQuantity; i++) {
                        if (i == bits[0]) {
                            currOp = tensorProd(currOp, this.operators[7]); //T
                        } else {
                            currOp = tensorProd(currOp, this.operators[1]); //I
                        }
                    }
                    
                    currOp.Multiply(this.stateVector, this.stateVector);
                    break;
                case opCodes.X:
                    //Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});

                    for (int i = 0; i < this.qbitQuantity; i++) {
                        if (i == bits[0]) {
                            currOp = tensorProd(currOp, this.operators[8]); //X
                        } else {
                            currOp = tensorProd(currOp, this.operators[1]); //I
                        }
                    }
                    
                    currOp.Multiply(this.stateVector, this.stateVector);
                    break;
                case opCodes.Y:
                    //Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});

                    for (int i = 0; i < this.qbitQuantity; i++) {
                        if (i == bits[0]) {
                            currOp = tensorProd(currOp, this.operators[9]); //Y
                        } else {
                            currOp = tensorProd(currOp, this.operators[1]); //I
                        }
                    }
                    
                    currOp.Multiply(this.stateVector, this.stateVector);
                    break;
                case opCodes.Z:
                    //Matrix<Complex> currOp = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});

                    for (int i = 0; i < this.qbitQuantity; i++) {
                        if (i == bits[0]) {
                            currOp = tensorProd(currOp, this.operators[10]); //Z
                        } else {
                            currOp = tensorProd(currOp, this.operators[1]); //I
                        }
                    }
                    
                    currOp.Multiply(this.stateVector, this.stateVector);
                    break;
                default:
                    break;
            }
        }

        public MathNet.Numerics.LinearAlgebra.Vector<Complex> getStateVector() {
            return this.stateVector;
        }

        public double[] getProbabilityVector() {
            double[] probVector = new double[(int) Math.Pow(2, this.qbitQuantity)];

            for (int i = 0; i < probVector.Length; i++) {
                probVector[i] = Math.Pow(Complex.Abs(this.stateVector[i]), 2);
            }

            return probVector;
        }

        //Private Functions//////////
        private Matrix<Complex> tensorProd(Matrix<Complex> matA, Matrix<Complex> matB) {
            Matrix<Complex> result = Matrix<Complex>.Build.Dense(matA.RowCount * matB.RowCount, matA.ColumnCount * matB.ColumnCount);

            for (int matARows = 0; matARows < matA.RowCount; matARows++) {
                for (int matACol = 0; matACol < matA.ColumnCount; matACol++) {
                    for (int matBRows = 0; matBRows < matB.RowCount; matBRows++) {
                        int resultsIndexRow = matARows * matB.RowCount + matBRows;
                        for (int matBCol = 0; matBCol < matB.ColumnCount; matBCol++) {
                            int resultsIndexCol = matACol * matB.ColumnCount + matBCol;

                            result[resultsIndexRow, resultsIndexCol] = matA[matARows, matACol] * matB[matBRows, matBCol];
                        }
                    }
                }
            }

            return result;
        }

        private MathNet.Numerics.LinearAlgebra.Vector<Complex> swapBits(MathNet.Numerics.LinearAlgebra.Vector<Complex> vec, int bitStart, int bitEnd) {
            MathNet.Numerics.LinearAlgebra.Vector<Complex> result = vec.Clone();

            int direction = (int) ((bitEnd - bitStart) / Math.Abs(bitEnd - bitStart));
            for (int i = 0; i < (int) Math.Abs(bitEnd - bitStart); i++) {
                Matrix<Complex> swapMat = genSwapMat(bitStart + (direction * i), vec.Count, (direction < 0));

                result = swapMat.Multiply(result);
            }

            return result;
        }

        //Generates the swap matrix, swapping the qbit with the adjacent qbit
        private Matrix<Complex> genSwapMat(int targetBit, int size, bool reverse = false) {
            Matrix<Complex> swapMat = Matrix<Complex>.Build.DenseOfArray(new Complex[,] {{1}});

            if (reverse) {
                for (int i = 0; i < size; i++) {
                    if (i == (targetBit - 1)) {
                        swapMat = tensorProd(swapMat, this.operators[6]);   //Swap
                    } else {
                        swapMat = tensorProd(swapMat, this.operators[1]);   //I
                    }
                }
            } else {
                for (int i = 0; i < size; i++) {
                    if (i == targetBit) {
                        swapMat = tensorProd(swapMat, this.operators[6]);   //Swap
                    } else {
                        swapMat = tensorProd(swapMat, this.operators[1]);   //I
                    }
                }
            }

            return swapMat;
        }

        private MathNet.Numerics.LinearAlgebra.Vector<Complex> measureBit(MathNet.Numerics.LinearAlgebra.Vector<Complex> vec, int bitNum) {
            MathNet.Numerics.LinearAlgebra.Vector<Complex> result;

            double probTrue = calcProbOfState(vec, bitNum);
            double probFalse = 1 - probTrue;

            double randomNum = this.rand.NextDouble();
            bool measuredState = (randomNum < probTrue);

            result = setNewMeasuredProb(vec, bitNum, measuredState ? probTrue : probFalse, measuredState);

            return result;
        }

        private double calcProbOfState(MathNet.Numerics.LinearAlgebra.Vector<Complex> vec, int bitNum, bool state = true) {
            double prob = 0;
            int bitMask = 1 << bitNum;

            for (int i = 0; i < vec.Count; i++) {
                if (state && (i & bitMask) != 0) {
                    prob += Math.Pow(Complex.Abs(vec[i]), 2);
                } else if (!state && (i & bitMask) == 0) {
                    prob += Math.Pow(Complex.Abs(vec[i]), 2);
                }
            }

            return prob;
        }

        private MathNet.Numerics.LinearAlgebra.Vector<Complex> setNewMeasuredProb(MathNet.Numerics.LinearAlgebra.Vector<Complex> vec, int bitNum, double probability, bool state = true) {
            MathNet.Numerics.LinearAlgebra.Vector<Complex> result = vec.Clone();
            int bitMask = 1 << bitNum;

            for (int i = 0; i < vec.Count; i++) {
                if (state && (i & bitMask) != 0) {
                    result[i] = (1 / Math.Sqrt(probability)) * vec[i];
                } else if (!state && (i & bitMask) == 0) {
                    result[i] = (1 / Math.Sqrt(probability)) * vec[i];
                } else {
                    result[i] = 0;
                }
            }

            return result;
        }
    }

}






