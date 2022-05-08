n = 6
m = 6
arr = [[0] * m for i in range(n)] #리스트 내포
a = n * m
i = 0
cnt = 0
while(1):
    for x in range(i, m -1 - i): #ㄱㄹ
        arr[i][x] =  a- cnt
        cnt += 1
    for y in range(i, n - 1 - i): #ㅅㄹ
        arr[y][m - 1 - i] = a - cnt #세로지만 여기선 열을 집어줘야함
        cnt += 1
    for z in range(m - 1 - i, i, -1):
        arr[n - 1 - i][z] = a - cnt #가로지만 여기선 행을 집어줘야함
        cnt += 1
    for w in range(n - 1 - i, i, -1):
        arr[w][i] = a - cnt
        cnt += 1

    i += 1
    if i == (m * n) // 2:
        break

for i in arr:
    print(i)

        