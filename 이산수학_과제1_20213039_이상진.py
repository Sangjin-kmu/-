# 입력
n = int(input("정방행렬의 차수를 입력하세요 : "))
mat = []
for i in range(n):
    row = list(map(int, input(f"{i+1}행 : ").split()))
    mat.append(row)

# 행렬식 (재귀)
def determinant(m):
    n = len(m)
    if n == 1: return m[0][0]
    if n == 2: return m[0][0]*m[1][1] - m[0][1]*m[1][0]
    d = 0
    for j in range(n):
        sub = [row[:j]+row[j+1:] for row in m[1:]]
        d += ((-1)**j) * m[0][j] * determinant(sub)
    return d

# 행렬식 이용한 역행렬
def inverse_by_determinant(m):
    n = len(m)
    d = determinant(m)
    if d == 0: return None
    cof = []
    for i in range(n):
        cof_row = []
        for j in range(n):
            sub = [row[:j]+row[j+1:] for k,row in enumerate(m) if k != i]
            cof_row.append(((-1)**(i+j)) * determinant(sub))
        cof.append(cof_row)
    adj = [[cof[j][i] for j in range(n)] for i in range(n)]  # 전치
    inv = [[adj[i][j]/d for j in range(n)] for i in range(n)]
    return inv

# 가우스-조던
def inverse_by_gauss_jordan(m):
    n = len(m)
    a = [row[:] for row in m]
    idn = [[float(i==j) for j in range(n)] for i in range(n)]
    for i in range(n):
        if a[i][i] == 0:  # 피벗 0이면 교환
            sw = None
            for k in range(i+1, n):
                if a[k][i] != 0:
                    sw = k
                    break
            if sw is None: return None
            a[i], a[sw] = a[sw], a[i]
            idn[i], idn[sw] = idn[sw], idn[i]
        piv = a[i][i]
        if piv == 0: return None
        for j in range(n):  # 피벗을 1로
            a[i][j] /= piv
            idn[i][j] /= piv
        for k in range(n):  # 다른 행 소거
            if k == i: continue
            ratio = a[k][i]
            for j in range(n):
                a[k][j] -= ratio * a[i][j]
                idn[k][j] -= ratio * idn[i][j]
    return idn

# 행렬 비교
def compare_matrices(m1, m2):
    eps=0.0000000001
    n = len(m1)
    for i in range(n):
        for j in range(n):
            if abs(m1[i][j]-m2[i][j]) > eps:
                return False
    return True

# 행렬 곱셈
def multiply_matrix(A, B):
    n = len(A)
    m = len(B[0])
    res = [[0]*m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = 0
            for k in range(len(B)):
                s += A[i][k] * B[k][j]
            res[i][j] = s
    return res

# 출력
r1 = inverse_by_determinant(mat)
if r1 is None:
    print("행렬식으로는 역행렬이 존재하지 않습니다.")
else:
    print("\n행렬식으로 구한 역행렬 :")
    for r in r1:
        print(" ".join(f"{v:8.4f}" for v in r))

r2 = inverse_by_gauss_jordan(mat)
if r2 is None:
    print("가우스-조던으로는 역행렬이 존재하지 않습니다.")
else:
    print("\n가우스-조던 소거법으로 구한 역행렬 :")
    for r in r2:
        print(" ".join(f"{v:8.4f}" for v in r))

if r1 and r2:
    if compare_matrices(r1, r2):
        print("\n두 방법의 결과가 동일합니다.")
    else:
        print("\n두 방법의 결과가 다릅니다.")

# 역행렬 검증 (A * A^-1 = I ?)
if r1:
    check = multiply_matrix(mat, r1)
    print("\n원래 행렬 × 역행렬 (행렬식 방식) :")
    for r in check:
        print(" ".join(f"{v:8.4f}" for v in r))
