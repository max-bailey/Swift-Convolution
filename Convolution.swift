import Foundation
import Accelerate
import Darwin
import Dispatch

/* Frequency-domain convolution */
func convolve(inputAudio: [Float], impulse: [Float]) -> [Float] {
    
    /* Convolve the two discrete-time signals */
    // Create an array containing the signal with extra space for the convolution to proagate
    // Pad zeros of length of other signal minus one
    var x: [Float] = inputAudio
    var h: [Float] = impulse
    
    // Find next power of two after signal lengths (convolution requires signal lengths of powers of two)
    let maxSignalLength: Int = inputAudio.count + impulse.count
    let nextPowerOfTwo: Int = Int(pow(2, ceil(log(Double(maxSignalLength)) / log(2))))
    // Append signal lengths to be next power of two
    let toAppendOnX: Int = nextPowerOfTwo - x.count
    let toAppendOnH: Int = nextPowerOfTwo - h.count
    let xAppendZeros: [Float] = [Float](repeating: 0, count: toAppendOnX)
    let hAppendZeros: [Float] = [Float](repeating: 0, count: toAppendOnH)
    x.append(contentsOf: xAppendZeros)
    h.append(contentsOf: hAppendZeros)
    
    // Perform FFT of 'inputAudio' and 'impulse'
    let fftInputAudio = fft(input: x)
    let fftImpulse = fft(input: h)
    
    // Element-wise multiply the two signals
    // Real
    var multipliedReal = [Float](repeating: 0, count: fftInputAudio.0.count)
    
    if fftInputAudio.0.count == fftImpulse.0.count {
        vDSP_vmul(fftInputAudio.0, 1, fftImpulse.0, 1, &multipliedReal, 1, vDSP_Length(fftImpulse.0.count))
    } else {
        fatalError("FFT of REAL input audio and impulse response are not equal.\n")
    }
    
    // Imaginary
    var multipliedImag = [Float](repeating: 0, count: fftInputAudio.1.count)
    
    if fftImpulse.1.count == fftImpulse.1.count {
        vDSP_vmul(fftInputAudio.1, 1, fftImpulse.1, 1, &multipliedImag, 1, vDSP_Length(fftImpulse.1.count))
    } else {
        fatalError("FFT of IMAGINARY input audio and impulse response are not equal.\n")
    }
    
    // Perform IFFT (INVERSE Fast Fourier Transform) of multipliedFFT
    let convolvedOutput: [Float] = ifft(real: multipliedReal, imaginary: multipliedImag)
    
    return convolvedOutput
}

/* FFT */
// From https://gist.github.com/jeremycochoy/45346cbfe507ee9cb96a08c049dfd34f
func fft(input: [Float]) -> ([Float], [Float]) {
    
    // The length of the input
    let length = vDSP_Length(input.count)
    
    // The power of two of two times the length of the input
    let log2n = vDSP_Length(ceil(log2(Float(length * 2))))
    
    // Create the instance of the FFT class which allows computing FFT of complex vector with length up to 'length'
    let fftSetup = vDSP.FFT(log2n: log2n, radix: .radix2, ofType: DSPSplitComplex.self)!
    
    // Input + output arrays
    var forwardInputReal = [Float](input) // Copy the signal
    var forwardInputImag = [Float](repeating: 0, count: Int(length))
    var forwardOutputReal = [Float](repeating: 0, count: Int(length))
    var forwardOutputImag = [Float](repeating: 0, count: Int(length))
    //var magnitudes = [Float](repeating: 0, count: Int(length)) // Magnitude of each frequency bin
    
    // Compute FFT
    forwardInputReal.withUnsafeMutableBufferPointer { forwardInputRealPtr in
        forwardInputImag.withUnsafeMutableBufferPointer { forwardInputImagPtr in
            forwardOutputReal.withUnsafeMutableBufferPointer { forwardOutputRealPtr in
                forwardOutputImag.withUnsafeMutableBufferPointer { forwardOutputImagPtr in
                    // Input
                    let forwardInput = DSPSplitComplex(realp: forwardInputRealPtr.baseAddress!, imagp: forwardInputImagPtr.baseAddress!)
                    // Output
                    var forwardOutput = DSPSplitComplex(realp: forwardOutputRealPtr.baseAddress!, imagp: forwardOutputImagPtr.baseAddress!)
                    
                    fftSetup.forward(input: forwardInput, output: &forwardOutput) // Call to the fft
                    //vDSP.absolute(forwardOutput, result: &magnitudes) // Compute the magitude of each frequency bin
                }
            }
        }
    }
    
    return (forwardOutputReal, forwardOutputImag)
}

func ifft(real: [Float], imaginary: [Float]) -> [Float] {
    
    // The length of the input
    var length: vDSP_Length!
    if real.count == imaginary.count {
        length = vDSP_Length(real.count)
    } else {
        fatalError("Real and imaginary signals are not equal length.\n")
    }
    
    // The power of two of two times the length of the input
    let log2n = vDSP_Length(ceil(log2(Float(length * 2))))
    
    // Create the instance of the FFT class which allows computing FFT of complex vector with length up to 'length'
    let fftSetup = vDSP.FFT(log2n: log2n, radix: .radix2, ofType: DSPSplitComplex.self)!
    
    // Input + output arrays
    var forwardInputReal = [Float](real) // Copy the signal
    var forwardInputImag = [Float](imaginary)
    var forwardOutputReal = [Float](repeating: 0, count: Int(length))
    var forwardOutputImag = [Float](repeating: 0, count: Int(length))
    var magnitudes = [Float](repeating: 0, count: Int(length)) // Magnitude of each frequency bin
    
    // Compute FFT
    forwardInputReal.withUnsafeMutableBufferPointer { forwardInputRealPtr in
        forwardInputImag.withUnsafeMutableBufferPointer { forwardInputImagPtr in
            forwardOutputReal.withUnsafeMutableBufferPointer { forwardOutputRealPtr in
                forwardOutputImag.withUnsafeMutableBufferPointer { forwardOutputImagPtr in
                    // Input
                    let forwardInput = DSPSplitComplex(realp: forwardInputRealPtr.baseAddress!, imagp: forwardInputImagPtr.baseAddress!)
                    // Output
                    var forwardOutput = DSPSplitComplex(realp: forwardOutputRealPtr.baseAddress!, imagp: forwardOutputImagPtr.baseAddress!)
                    
                    fftSetup.inverse(input: forwardInput, output: &forwardOutput) // Call to the fft
                    vDSP.absolute(forwardOutput, result: &magnitudes) // Compute the magitude of each frequency bin
                }
            }
        }
    }
    
    return forwardOutputReal
}

