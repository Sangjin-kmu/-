def print_matrix(matrix, title=None):
    if title:
        print(f"\n{title}")
    for row in matrix:
        print(row)
    print()

def input_relation_matrix():
    print("관계행렬의 크기를 입력하세요:")
    n = int(input("n = "))
    print("각 행을 공백으로 구분해 입력하세요 (0 또는 1):")
    matrix = []
    for i in range(n):
        row = list(map(int, input(f"{i+1}행: ").split()))
        if len(row) != n:
            raise ValueError("입력한 원소의 개수가 올바르지 않습니다.")
        matrix.append(row)
    A = list(range(1, n+1))
    return matrix, A

def is_reflexive(R):
    n = len(R)
    return all(R[i][i] == 1 for i in range(n))


def is_symmetric(R):
    n = len(R)
    return all(R[i][j] == R[j][i] for i in range(n) for j in range(n))


def is_transitive(R):
    n = len(R)
    for i in range(n):
        for j in range(n):
            if R[i][j]:
                for k in range(n):
                    if R[j][k] and not R[i][k]:
                        return False
    return True


def is_equivalence(R):
    return is_reflexive(R) and is_symmetric(R) and is_transitive(R)

def reflexive_closure(R):
    n = len(R)
    newR = [row[:] for row in R]
    for i in range(n):
        newR[i][i] = 1
    return newR


def symmetric_closure(R):
    n = len(R)
    newR = [row[:] for row in R]
    for i in range(n):
        for j in range(n):
            if R[i][j] == 1:
                newR[j][i] = 1
    return newR


def transitive_closure(R, show_steps=False):
    n = len(R)
    closure = [row[:] for row in R]
    changed = True
    step = 1
    while changed:
        changed = False
        added = []

        for i in range(n):
            for j in range(n):
                if closure[i][j] == 0:
                    for k in range(n):
                        if closure[i][k] and closure[k][j]:
                            closure[i][j] = 1
                            changed = True
                            if show_steps:
                                added.append((i+1, k+1, j+1))
                            break

        if show_steps and added:
            print(f"\n[추이 확장 step {step}] 추가된 경로:")
            for a, b, c in added:
                print(f"{a} → {b} → {c} 경로로 인해 {a} → {c} 추가됨")
        step += 1

    if show_steps:
        print_matrix(closure, "최종 연결관계 행렬:")

    return closure


def warshall_closure(R):
    n = len(R)
    W = [row[:] for row in R]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                W[i][j] = int(W[i][j] or (W[i][k] and W[k][j]))
    return W


def show_property_results(R):
    ref = is_reflexive(R)
    sym = is_symmetric(R)
    tra = is_transitive(R)
    print("\n[관계 성질 판별 결과]")
    print(f"  반사성 (Reflexive): {'true' if ref else 'false'}")
    print(f"  대칭성 (Symmetric): {'true' if sym else 'false'}")
    print(f"  추이성 (Transitive): {'true' if tra else 'false'}")
    print("")
    return ref, sym, tra


def print_equivalence_classes(R, A):
    print("\n[동치류 출력]")
    for i in range(len(A)):
        cls = [A[j] for j in range(len(A)) if R[i][j] == 1]
        print(f"[{A[i]}] = {cls}")


def main():
    R, A = input_relation_matrix()
    print_matrix(R, "입력된 관계행렬")

    show_property_results(R)
    if is_equivalence(R):
        print("\n이 관계는 동치 관계입니다.")
        print_equivalence_classes(R, A)

    print("반사 폐포 (Reflexive Closure)")
    print_matrix(R, "변환 전:")
    R_ref = reflexive_closure(R)
    print_matrix(R_ref, "변환 후:")
    show_property_results(R_ref)
    if is_equivalence(R_ref):
        print("반사 폐포 후 이 관계는 동치 관계입니다.")
        print_equivalence_classes(R_ref, A)

    print("대칭 폐포 (Symmetric Closure)")
    print_matrix(R, "변환 전:")
    R_sym = symmetric_closure(R)
    print_matrix(R_sym, "변환 후:")
    show_property_results(R_sym)
    if is_equivalence(R_sym):
        print("대칭 폐포 후 이 관계는 동치 관계입니다.")
        print_equivalence_classes(R_sym, A)

    print("추이 폐포 (Transitive Closure)")
    print_matrix(R, "변환 전:")
    R_tra = transitive_closure(R)
    show_property_results(R_tra)
    if is_equivalence(R_tra):
        print("추이 폐포 후 이 관계는 동치 관계입니다.")
        print_equivalence_classes(R_tra, A)

    print("-모든 폐포 순차 적용 (Reflexive → Symmetric → Transitive)")
    step1 = reflexive_closure(R)
    print_matrix(step1, "-반사 폐포 적용 후:")
    show_property_results(step1)
    if is_equivalence(step1):
        print("반사 폐포 후 동치 관계입니다.")
        print_equivalence_classes(step1, A)

    step2 = symmetric_closure(step1)
    print_matrix(step2, "-대칭 폐포 적용 후:")
    show_property_results(step2)
    if is_equivalence(step2):
        print("대칭 폐포 후 동치 관계입니다.")
        print_equivalence_classes(step2, A)

    step3 = transitive_closure(step2)
    print_matrix(step3, "-추이 폐포 적용 후:")
    show_property_results(step3)
    if is_equivalence(step3):
        print("모든 폐포 적용 후 동치 관계입니다.")
        print_equivalence_classes(step3, A)
    else:
        print("모든 폐포를 적용해도 동치 관계가 아닙니다.")

    print("[추가구현] Warshall 알고리즘으로 계산한 추이 폐포")
    W = warshall_closure(R)
    print_matrix(W, "Warshall 결과:")


if __name__ == "__main__":
    main()
