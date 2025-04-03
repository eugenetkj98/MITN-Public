using Distributed
workers()
rmprocs(workers())
addprocs(3)
@distributed for i in 1:10
    Threads.@threads for j in 1:30
        println("Worker ID: $(Distributed.myid()), on thread: $(Threads.threadid())")
        sleep(1)
    end
end