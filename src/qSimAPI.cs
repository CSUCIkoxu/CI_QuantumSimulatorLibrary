using System.Collections.Generic;
using System;
using System.IO;
using System.IO.Pipes;
using System.Diagnostics;

namespace qSim {
    public enum opCodes {
        nop = 0,
        M = 1,
        H = 2,
        I = 3,
        CNOT = 4,
        CCNOT = 5,
        R = 6,
        S = 7,
        SWAP = 8,
        T = 9,
        X = 10,
        Y = 11,
        Z = 12
    }

    public enum PauliAxis {
        I = 0,
        X = 1,
        Y = 2,
        Z = 3
    }
    
    public class qSimAPI {
        private qSimulator qSim;
        private int bitNum;
        private List<opCodes> opList;
        private List<int[]> qbitList;
        private List<double[]> paramList;

        public qSimAPI() {
            this.qSim = new qSimulator();
            this.bitNum = 0;
            this.opList = new List<opCodes>{};
            this.qbitList = new List<int[]>{};
            this.paramList = new List<double[]>{};
        }

        public qSimAPI(int startingBits) {
            this.qSim = new qSimulator(startingBits);
            this.bitNum = startingBits;
            this.opList = new List<opCodes>{};
            this.qbitList = new List<int[]>{};
            this.paramList = new List<double[]>{};
        }

        public qSimAPI(int startingBits, int seed) {
            this.qSim = new qSimulator(startingBits, seed);
            this.bitNum = startingBits;
            this.opList = new List<opCodes>{};
            this.qbitList = new List<int[]>{};
            this.paramList = new List<double[]>{};
        }


        public void setBitNum(int newBitNum) {bitNum = newBitNum;this.qSim.newShape(newBitNum);}
        public int getBitNum() {return bitNum;}

        //DEBUG ONLY!!!
        public void setInputBuffer(opCodes[] newOpCodes, int[][] newQBits, double[][] newParams) {
            this.opList = new List<opCodes>(newOpCodes);
            this.qbitList = new List<int[]> {};
            this.paramList = new List<double[]> {};

            for (int i = 0; i < newQBits.Length; i++) {
                this.qbitList[i] = newQBits[i];
            }
            for (int i = 0; i < newParams.Length; i++) {
                this.paramList[i] = newParams[i];
            }
        }

        public void pushOp(opCodes op, int[] qBits, double[] opParams) {
            this.opList.Add(op);
            this.qbitList.Add(qBits);
            this.paramList.Add(opParams);
        }

        public (opCodes, int[], double[]) popOp() {
            int lastIndex = opList.Count-1;

            opCodes op = this.opList[lastIndex];
            int[] bits = this.qbitList[lastIndex];
            double[] parameters = this.paramList[lastIndex];

            (opCodes, int[], double[]) popped = (op, bits, parameters);
            this.opList.RemoveAt(lastIndex);
            this.qbitList.RemoveAt(lastIndex);
            this.paramList.RemoveAt(lastIndex);

            //(opCodes, long, long, long, double) poppedOut = ((opCodes) popped.Item1, popped.Item2, popped.Item3, popped.Item4, popped.Item5);

            return popped;
        }

        public (opCodes, int[], double[]) removeOp(int index) {
            opCodes op = this.opList[index];
            int[] bits = this.qbitList[index];
            double[] parameters = this.paramList[index];

            (opCodes, int[], double[]) popped = (op, bits, parameters);
            this.opList.RemoveAt(index);
            this.qbitList.RemoveAt(index);
            this.paramList.RemoveAt(index);

            //(opCodes, long, long, long, double) poppedOut = ((opCodes) popped.Item1, popped.Item2, popped.Item3, popped.Item4, popped.Item5);

            return popped;
        }

        public void insertOp(int index, opCodes op, int[] qBits, double[] opParams) {
            this.opList.Insert(index, op);
            this.qbitList.Insert(index, qBits);
            this.paramList.Insert(index, opParams);
        }

        public void pushM(int bitIndex) {
            pushOp(opCodes.M, new int[] {bitIndex}, new double[0]);
        }
        public void pushH(int bitIndex) {
            pushOp(opCodes.H, new int[] {bitIndex}, new double[0]);
        }
        public void pushI(int bitIndex) {
            pushOp(opCodes.I, new int[] {bitIndex}, new double[0]);
        }
        public void pushCNOT(int control, int target) {
            pushOp(opCodes.CNOT, new int[] {control, target}, new double[0]);
        }
        public void pushCCNOT(int control1, int control2, int target) {
            pushOp(opCodes.CCNOT, new int[] {control1, control2, target}, new double[0]);
        }
        public void pushR(int bitIndex, PauliAxis pauliAxis, double theta) {
            pushOp(opCodes.R, new int[] {bitIndex}, new double[] {(double) pauliAxis, theta});
        }
        public void pushS(int bitIndex) {
            pushOp(opCodes.S, new int[] {bitIndex}, new double[0]);
        }
        public void pushSWAP(int bitIndex, int swapIndex) {
            pushOp(opCodes.SWAP, new int[] {bitIndex, swapIndex}, new double[0]);
        }
        public void pushT(int bitIndex) {
            pushOp(opCodes.T, new int[] {bitIndex}, new double[0]);
        }
        public void pushX(int bitIndex) {
            pushOp(opCodes.X, new int[] {bitIndex}, new double[0]);
        }
        public void pushY(int bitIndex) {
            pushOp(opCodes.Y, new int[] {bitIndex}, new double[0]);
        }
        public void pushZ(int bitIndex) {
            pushOp(opCodes.Z, new int[] {bitIndex}, new double[0]);
        }

        public void pushHOnAll() {
            for (int i = 0; i < bitNum; i++) {
                pushH(i);
            }
        }

        public double[] run() {
            double[] results;

            for (int i = 0; i < opList.Count; i++) {
                this.qSim.runOp(this.opList[i], this.qbitList[i], this.paramList[i]);
            }

            results = this.qSim.getProbabilityVector();

            return results;
        }
    }
}