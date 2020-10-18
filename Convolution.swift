import Foundation
import Accelerate

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
