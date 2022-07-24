
Public Function FFT(ByVal A As Complex()) As Complex()
    Dim n As Integer = A.Length
    If n = 1 Then
        Return A
    Else
        Dim B(n - 1) As Complex
        For k As Integer = 0 To n \ 2 - 1
            B(k) = A(k) + A(n \ 2 + k)
            B(n \ 2 + k) = A(k) - A(n \ 2 + k)
        Next
        Dim C As Complex() = FFT(B(0 To n \ 2 - 1))
        Dim D As Complex() = FFT(B(n \ 2 To n - 1))
        For k As Integer = 0 To n \ 2 - 1
            B(k) = C(k) + Complex.Exp(-2 * Math.PI * k * 1 / n) * D(k)
            B(n \ 2 + k) = C(k) - Complex.Exp(-2 * Math.PI * k * 1 / n) * D(k)
        Next
        Return B
    End If
End Function
