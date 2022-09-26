

FFT(A())

n = UBound(A)

If n = 1 Then
    B(0) = A(0)
Else
    For k = 0 To n \ 2 - 1
        B(k) = A(k) + A(n \ 2 + k)
        B(n \ 2 + k) = A(k) - A(n \ 2 + k)
    Next k
    
    C = FFT(B(0 To n \ 2 - 1))
    D = FFT(B(n \ 2 To n - 1))
    
    For k = 0 To n \ 2 - 1
        B(k) = C(k) + Exp(-2 * pi * k * i / n) * D(k)
        B(n \ 2 + k) = C(k) - Exp(-2 * pi * k * i / n) * D(k)
    Next k
End If

FFT = B
End FFT
