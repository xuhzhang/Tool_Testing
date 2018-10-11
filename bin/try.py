from concurrent.futures import ProcessPoolExecutor
import os, time, random

def task(n, m):
    print("%s is running" % os.getpid())
    time.sleep(2)
    return(str(n**2) + "," + m + ",")

if __name__ == "__main__":
    #p = ProcessPoolExecutor(4)
    p = ProcessPoolExecutor(10)
    #l = []
    start=time.time()
    obj = ''
    for i in range(10):
        obj += p.submit(task, i, "a").result()
        #l.append(obj)

    p.shutdown()
    print('='*30)
    print(obj)
    print(time.time() - start)
