f(x)= a*x**b
#f(x) = 1/(A2 * x**2 + A4 * x**4)
a = 100
b = -1
#A2=1
#A4=1
fit [1:10] f(x) "ft_L101_W21.txt" u 1:3:4 via a,b
ti = sprintf("%.2fx**%.2f", a, b)
set logscale x
set logscale y
plot "ft_L101_W21.txt" u 1:3:4 with errorlines lt rgb "#ff0000" title "L101 W21", \
f(x) with lines lt rgb "#ff00ff" title ti
