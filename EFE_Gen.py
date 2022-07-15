import sympy as sp
import numpy as np
import sys
import subprocess, os
import tex

# sp.init_printing()

def get_christoffel_symbols(metric, axes):
    metric_inv = metric.inv()
    # print('Inverse Metric')
    # print(metric_inv)
    # christoffel = sp.Matrix.zeros(4, 4, 4)
    christoffel = np.zeros([4, 4, 4], dtype = type(sp.Symbol('')))

    for i in range(4):
        for j in range(4):
            for k in range(4):
                # print('Chrostoffel', i, j, k)
                for s in range(4):
                    # print(metric_inv[s, i] * (sp.diff(metric[s, j], axes[k]) + sp.diff(metric[s, k], axes[j]) - sp.diff(metric[j, k], axes[i])))
                    christoffel[i][j][k] += metric_inv[s, i] * (sp.diff(metric[s, j], axes[k]) + sp.diff(metric[s, k], axes[j]) - sp.diff(metric[j, k], axes[i]))
                christoffel[i][j][k] = christoffel[i][j][k] / 2
                christoffel[i][j][k] = sp.simplify(christoffel[i][j][k])
                    # Chrostoffel[i][j][k] += 
                # print(christoffel[i][j][k])
    # print('Christoffel Symbols: \n', christoffel)
    return christoffel

def get_reimann_tensor(christoffel_symbols, axes):
    reimann = np.zeros([4, 4, 4, 4], dtype = type(sp.Symbol('')))
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    differential_part = 0
                    coeff_sum_part = 0
                    differential_part = sp.diff(christoffel_symbols[i][j][l], axes[k]) - sp.diff(christoffel_symbols[i][k][l], axes[j])
                    for p in range(4):
                        coeff_sum_part += christoffel_symbols[p][j][l] * christoffel_symbols[i][p][k] - christoffel_symbols[p][k][l] * christoffel_symbols[i][p][j]
                    reimann[i][j][k][l] = differential_part + coeff_sum_part
                    reimann[i][j][k][l] = sp.simplify(reimann[i][j][k][l])
    
    # reimann = sp.simplify(reimann)
    # print('\nReimann Curvature Tensor: \n', reimann)
    return reimann

def get_ricci_tensor(reimann_tensor):
    ricci = np.zeros([4, 4], dtype=type(sp.Symbol('')))
    for i in range(4):
        for j in range(4):
            for k in range(4):
                ricci[i][j] += reimann_tensor[k, i, k, j]
            ricci[i][j] = sp.simplify(ricci[i][j])

    # ricci = sp.simplify(ricci)
    print('\nRicci Curvature Tensor: \n', sp.Matrix(ricci))
    return sp.Matrix(ricci)

def raise_one_index(tensor, metric):
    shape = tensor.shape
    rank = len(shape)
    metric_inv = metric.inv()
    raised_tensor = np.zeros(shape, dtype=type(sp.Symbol('')))
    for i in range(4):
        for j in range(4):
            for k in range(4):
                raised_tensor[i][j] += metric_inv[i, k] * tensor[k, j]
            raised_tensor[i][j] = sp.simplify(raised_tensor[i][j])
    # raised_tensor = sp.simplify(raised_tensor)
    # print('\n Raised Tensor: \n', raised_tensor)
    return sp.Matrix(raised_tensor)
    

def get_curvature_scalar(raised_ricci):
    curvature_scalar = 0
    for i in range(4):
        for j in range(4):
            curvature_scalar += raised_ricci[i, j]
    curvature_scalar = sp.simplify(curvature_scalar)
    print('\nCurvature Scalar: \n', curvature_scalar)
    return curvature_scalar

def conform_compacted_metric(axes):
    return sp.Matrix(([-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, sp.sin(axes[1])**2, 0], [0, 0, 0, sp.sin(axes[1])**2 * sp.sin(axes[2])**2]))

def spherical_metric(axes):
    c = sp.Symbol('c')
    return sp.Matrix(([-c**2, 0, 0, 0], [0, 1, 0, 0], [0, 0, axes[1]**2, 0], [0, 0, 0, axes[1]**2 * sp.sin(axes[2]) ** 2]))

def FRW_metric(axes):
    a, k = sp.symbols('a k')
    return sp.Matrix(([-1, 0, 0, 0],[0, a**2/(1-k*axes[1]**2), 0, 0], [0, 0, a**2*axes[1]**2, 0],[0, 0, 0, a**2*axes[1]**2*sp.sin(axes[2])**2]))

# def write_into_latex_file():


def main():
    x = sp.symbols('x0 x1 x2 x3')
    sp.init_printing(use_unicode=True)

    c = sp.Symbol('c')

    metric_tensor = conform_compacted_metric(x)

    christoffel_symbols = get_christoffel_symbols(metric_tensor, x)
    reimann_curvature_tensor = get_reimann_tensor(christoffel_symbols, x)
    ricci_curvature_tensor = get_ricci_tensor(reimann_curvature_tensor)
    raised_ricci_tensor = raise_one_index(ricci_curvature_tensor, metric_tensor)
    curvature_scalar = get_curvature_scalar(raised_ricci_tensor)

    einstein_tensor = ricci_curvature_tensor - curvature_scalar * metric_tensor
    print('\n\nEinstein Tensor: \n\n', einstein_tensor)

    filename = "EinsteinFieldEquations"
    fileHandle = open(filename + ".tex", 'w')
    with fileHandle as file:
        file.write('\\documentclass{article}\n')
        file.write('\\usepackage{amsmath}\n')
        # file.write('\\usepackage[utf8]{inputenc}\n')
        file.write('\\title{General Relativity Assignment}\n')
        file.write('\\author{Nimesh Khandelwal}\n')
        file.write('\\begin{document}\n')
        file.write('\\maketitle\n')
        file.write('\\section{Important Symbols and Tensors}\n')
        file.write('\\subsection{Metric}\n')
        file.write('$$ g_\\mu{_\\nu} = ' + sp.latex(metric_tensor) + '$$\n')
        file.write('\\subsection{Christoffel Symbols}\n')
        # file.write('$$' + sp.latex(christoffel_symbols) + '$$\n')
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    if christoffel_symbols[i][j][k] == 0:
                        continue
                    file.write('$$ \Gamma^' + str(i) + '{_' + str(j) + '{_' + str(k) + '}} = ' + sp.latex(christoffel_symbols[i][j][k]) + '$$\n')
        file.write('\\subsection{Ricci Tensor}\n')
        file.write('$$ R_\\mu{_\\nu} = ' + sp.latex(ricci_curvature_tensor) + '$$\n')
        file.write('\\subsection{Curvature scalar}\n')
        file.write('$$ R = ' + sp.latex(curvature_scalar) + '$$\n')
        file.write('\\subsection{Einstein Tensor}\n')
        file.write('$$ G_\\mu{_\\nu} = R_\\mu{_\\nu} - Rg_\\mu{_\\nu} = ' + sp.latex(einstein_tensor) + '$$\n')
        file.write('\\end{document}\n')
    fileHandle.close()

    pdf_conversion_process = subprocess.call('pdflatex '+ filename, shell=True)

    if pdf_conversion_process != 0:
        print('Check result! Exit code not 0')
    else:
        os.system("xdg-open " + filename + ".pdf")



if __name__ == "__main__":
    main()
