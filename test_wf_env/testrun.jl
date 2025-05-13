println("This client is running with $(Threads.nthreads()) threads!")
N = 100


B = zeros(N,N)

for i in 1:N
    B[i,:] = rand(N)
end

println("I finished the calculation. The total of all entries is $(sum(B)) :D")