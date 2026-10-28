import os
import graphviz

# Ensure Graphviz binary directory (MSYS2 UCRT64 or Program Files) is on PATH
for gv_path in [r"C:\msys64\ucrt64\bin", r"C:\Program Files\Graphviz\bin"]:
    if os.path.isdir(gv_path) and gv_path not in os.environ.get("PATH", ""):
        os.environ["PATH"] = gv_path + os.pathsep + os.environ.get("PATH", "")


def build_global_flowchart():
    """Generates the overall program flowchart (`diagrama_graphviz.png`)."""
    dot = graphviz.Digraph(comment='Flowchart 2D Elastoplastic BEM (ELQUABE)', format='png')
    dot.attr(rankdir='LR', compound='true', dpi='300', nodesep='0.35', ranksep='0.45')

    node_attr = {
        'shape': 'box',
        'style': 'filled,rounded',
        'fillcolor': '#d4edda',
        'color': '#c3e6cb',
        'fontname': 'sans-serif',
        'fontsize': '11',
    }
    plast_attr = {
        'shape': 'box',
        'style': 'filled,rounded',
        'fillcolor': '#cce5ff',
        'color': '#b8daff',
        'fontname': 'sans-serif',
        'fontsize': '11',
    }
    circle_attr = {
        'shape': 'circle',
        'width': '0.15',
        'style': 'filled',
        'fillcolor': '#e8e8ff',
        'color': '#c3c6d8',
        'label': '',
    }
    oval_attr = {
        'shape': 'ellipse',
        'style': 'filled',
        'fillcolor': '#e8e8ff',
        'color': '#c3c6d8',
        'fontname': 'sans-serif',
        'fontsize': '11',
    }
    file_attr = {
        'shape': 'parallelogram',
        'style': 'filled',
        'fillcolor': '#fff3cd',
        'color': '#ffeeba',
        'fontname': 'sans-serif',
        'fontsize': '10',
    }

    dot.node('AA', 'BEGIN', **oval_attr)
    dot.node('G', 'END', **oval_attr)

    with dot.subgraph(name='cluster_main') as c:
        c.attr(
            label='MAIN PROGRAM (ELQUABE - 2D Elastoplastic BEM)',
            fontname='sans-serif',
            fontsize='12',
            style='dashed,rounded',
            color='gray',
            bgcolor='#fafafa',
        )

        c.node('BB', 'INPUT FILE\n(eq12cb.dat)', **file_attr)
        c.node('B', 'INPUTEQ\n(Mesh & Cells)', **node_attr)
        c.node('C', 'GHMATEQ\n(H, G, Ab)', **node_attr)
        c.node('D', 'SLNPD\n(Elastic Solver)', **node_attr)
        c.node('E', 'INTEREQ\n(Au, Asig, Stress)', **node_attr)

        c.node('H', 'NL_COEFS\n(2D Cells: Eb, Eu, Esig)', **plast_attr)
        c.node('preINC', **circle_attr)
        c.node('I', 'NL_SOLVE / P_M_ITER\n(SIG_SCALE, PL_FLOW,\nSLNPD Δλ, UPDATEV)', **plast_attr)
        c.node('postINC', **circle_attr)

        c.node('F', 'OUTPTEQ\n(Write Results)', **node_attr)
        c.node('FF', 'OUTPUT FILES\n(eq12cb.out, dpbou/int.dat,\nstbou/int.dat)', **file_attr)

        with c.subgraph() as s_in:
            s_in.attr(rank='same')
            s_in.node('BB')
            s_in.node('B')

        with c.subgraph() as s_out:
            s_out.attr(rank='same')
            s_out.node('F')
            s_out.node('FF')

        c.edge('BB', 'B')
        c.edge('B', 'C')
        c.edge('C', 'D')
        c.edge('D', 'E')
        c.edge('E', 'H')
        c.edge('H', 'preINC')
        c.edge('preINC', 'I')
        c.edge('I', 'postINC')
        c.edge('postINC', 'F')
        c.edge('F', 'FF')

        c.edge(
            'postINC',
            'preINC',
            label='NEXT ITERATION /\nLOAD INCREMENT',
            fontname='sans-serif',
            fontsize='9',
        )

    dot.edge('AA', 'B', lhead='cluster_main')
    dot.edge('F', 'G', ltail='cluster_main')

    with dot.subgraph(name='cluster_legend') as leg:
        leg.attr(label='Legend', fontname='sans-serif', fontsize='11', color='black', style='solid', bgcolor='#ffffff')
        leg.node('leg_se', 'Start / End', **oval_attr)
        leg.node('leg_sub', 'Elastic / I/O Subroutine', **node_attr)
        leg.node('leg_pl', '2D Plasticity Subroutine', **plast_attr)
        leg.node('leg_file', 'Data / File', **file_attr)
        leg.edge('leg_se', 'leg_sub', style='invis')
        leg.edge('leg_sub', 'leg_pl', style='invis')
        leg.edge('leg_pl', 'leg_file', style='invis')

    dot.render('diagrama_graphviz', format='png', cleanup=True)
    print("O diagrama global foi exportado com sucesso como 'diagrama_graphviz.png'!")


def build_nonlinear_flowchart():
    """
    Generates a detailed flowchart (`diagrama_nl_graphviz.png`) of the 2D
    non-linear elastoplastic execution (`NL_COEFS` and `NL_SOLVE` / `P_M_ITER`).
    """
    dot = graphviz.Digraph(comment='Detailed 2D Nonlinear Elastoplastic BEM Flowchart', format='png')
    dot.attr(rankdir='TB', compound='true', dpi='300', nodesep='0.35', ranksep='0.42')

    oval_attr = {
        'shape': 'ellipse',
        'style': 'filled',
        'fillcolor': '#e8e8ff',
        'color': '#c3c6d8',
        'fontname': 'sans-serif',
        'fontsize': '10.5',
    }
    cell_attr = {
        'shape': 'box',
        'style': 'filled,rounded',
        'fillcolor': '#d4edda',
        'color': '#9fd3ab',
        'fontname': 'sans-serif',
        'fontsize': '10',
    }
    plast_attr = {
        'shape': 'box',
        'style': 'filled,rounded',
        'fillcolor': '#cce5ff',
        'color': '#99caff',
        'fontname': 'sans-serif',
        'fontsize': '10',
    }
    iter_attr = {
        'shape': 'box',
        'style': 'filled,rounded',
        'fillcolor': '#f8d7da',
        'color': '#f1b0b7',
        'fontname': 'sans-serif',
        'fontsize': '10',
    }
    dec_attr = {
        'shape': 'diamond',
        'style': 'filled',
        'fillcolor': '#fff3cd',
        'color': '#ebd386',
        'fontname': 'sans-serif',
        'fontsize': '9.5',
    }

    dot.node('START', 'From INTEREQ\n(Ab, Au, Asig, Elastic ESTRS)', **oval_attr)
    dot.node('END_NL', 'To OUTPTEQ\n(Final Displacements, Stresses & EPSTN)', **oval_attr)

    # 1. Subgraph: 2D Domain Cell Integration (NL_COEFS)
    with dot.subgraph(name='cluster_nl_coefs') as c1:
        c1.attr(
            label='1. 2D Domain Cell Integration — NL_COEFS (OpenMP Parallel over Nodes IP = 1..NTP)',
            fontname='sans-serif',
            fontsize='11.5',
            style='dashed,rounded',
            color='#28a745',
            bgcolor='#f7fcf8',
        )
        c1.node('INT_CELL', 'INT_CELL_2D\nLoop over 2D Cells ICELL = 1..NCELL\n(NODC = 6 Triangles or NODC = 8 Quads)', **cell_attr)
        c1.node('DEC_SING', 'Singular Node\nin Cell?\n(NCS > 0)', **dec_attr)
        c1.node('ADAPT', 'ADAPTINT_CELL_2D\nAdaptive Gauss Order (CHOSEGP_1D)\n& Sub-Element Subdivision', **cell_attr)
        c1.node('SING', 'SIN_CELL_2D & CELL_BOUND_2D\nPolar/Duffy Sub-Triangles +\nLine Integral (Eq. 9.10)', **cell_attr)
        c1.node('EVAL_KE', 'EVAL_KE_2D + Free Terms F_ijkl\nAssemble Initial Stress Matrices:\nEb_MAT (NTU=4), Eu_MAT (NTU=4), Esig_MAT (NTU=5)', **cell_attr)

        c1.edge('INT_CELL', 'DEC_SING')
        c1.edge('DEC_SING', 'ADAPT', label='No (Regular /\nNear-Singular)', fontname='sans-serif', fontsize='9')
        c1.edge('DEC_SING', 'SING', label='Yes (Singular\nVertex/Mid-node)', fontname='sans-serif', fontsize='9')
        c1.edge('ADAPT', 'EVAL_KE')
        c1.edge('SING', 'EVAL_KE')

    # 2. Subgraph: Coupled Plastic Influence Matrices & Initial Scaling (NL_SOLVE)
    with dot.subgraph(name='cluster_nl_init') as c2:
        c2.attr(
            label='2. Initial Yield Scaling & Coupled Influence Matrices — NL_SOLVE',
            fontname='sans-serif',
            fontsize='11.5',
            style='dashed,rounded',
            color='#0066cc',
            bgcolor='#f5faff',
        )
        c2.node('MAX_YIELD', 'PL_FLOW (All Nodes IP = 1..NTP)\nCompute Max Equivalent Stress EFSIG\n& Initial Yield Factor FSCAL = UNIAX / EFSIG', **plast_attr)
        c2.node('DEC_ELAST', 'EFSIG <= UNIAX?\n(Purely Elastic)', **dec_attr)
        c2.node('FORM_MAT', 'FORM_PLASTIC_MATRICES (In-Memory LAPACK DGESV)\n• AC_MAT  = AB_MAT⁻¹ · EB_MAT   (SLNPD_MULTI)\n• EUI_MAT = EU_MAT   + AU_MAT · AC_MAT\n• E_MAT   = ESIG_MAT + ASIG_MAT · AC_MAT', **plast_attr)
        c2.node('SCALE_INIT', 'Scale Elastic State to Initial Yield (FSCAL)\n& Set Load Increment Step WOVER (NINCS, FINCR)', **plast_attr)

        c2.edge('MAX_YIELD', 'DEC_ELAST')
        c2.edge('DEC_ELAST', 'SCALE_INIT', label='No (Plastic\nYielding)', fontname='sans-serif', fontsize='9')
        c2.edge('SCALE_INIT', 'FORM_MAT')

    # 3. Subgraph: Incremental-Iterative Plastic Multiplier Solver (P_M_ITER)
    with dot.subgraph(name='cluster_pm_iter') as c3:
        c3.attr(
            label='3. Incremental-Iterative Plastic Multiplier Solver — P_M_ITER',
            fontname='sans-serif',
            fontsize='11.5',
            style='dashed,rounded',
            color='#cc3344',
            bgcolor='#fff8f8',
        )
        c3.node('INC_START', 'Outer Load Increment Loop (INCRS = 1..NINCS)\nFSCAL = FSCAL + WOVER · FMULT', **plast_attr)
        c3.node('SIG_SCALE', 'SIG_SCALE & SIGCROSS\nApply Elastic Trial Increment, Find Yield Crossing ROVER,\nIdentify Active Yielded Nodes (NLAMD) & Excess Stress YE', **plast_attr)
        c3.node('PL_FLOW_INC', 'PL_FLOW & DF_DSIG\nEvaluate Yield Function & Flow Derivatives {∂f/∂σ}\n(Tresca, von Mises, Mohr-Coulomb, or Drucker-Prager)', **plast_attr)

        # Inner Iterative Loop
        c3.node('ITER_ASSM', 'Inner Iteration Loop (ITER = 1..MITER)\nDF_MATRX: Compute [De]{∂f/∂σ}, ψ, H_e\nAssemble Active System [A_λ]{Δλ} = {b_λ} (dim = NLAMD)', **iter_attr)
        c3.node('CONV_CHK', 'Convergence?\nRERR = ||r_λ|| / UNIAX\n< TOLER', **dec_attr)
        c3.node('SOLVE_UPD', 'SLNPD (LAPACK DGESV): Solve for Plastic Multipliers {δλ}\nUPDATEV: Update Boundary X_B (+AC_MAT·Δσᵖ),\nInternal Displacements DSOL (+EUI_MAT·Δσᵖ) & ESTRS (+E_MAT·Δσᵖ)', **iter_attr)

        c3.node('ACCUM_PL', 'Update Plastic State at Converged Increment\nAccumulate Equivalent Plastic Strain EPSTN (ε̄ᵖ)\n& Kinematic Back-Stress BACKS (α_ij)', **plast_attr)
        c3.node('MORE_INC', 'More Load\nIncrements?\n(INCRS < NINCS)', **dec_attr)
        c3.node('RECONSTR', 'Reconstruct Final Physical Fields\n• Boundary Displacements FI & Tractions DFI\n• Physical Stresses (ESTRS + BACKS) · 2G & von Mises', **plast_attr)

        c3.edge('INC_START', 'SIG_SCALE')
        c3.edge('SIG_SCALE', 'PL_FLOW_INC')
        c3.edge('PL_FLOW_INC', 'ITER_ASSM')
        c3.edge('ITER_ASSM', 'CONV_CHK')
        c3.edge('CONV_CHK', 'SOLVE_UPD', label='No (Next Iteration)', fontname='sans-serif', fontsize='9')
        c3.edge('SOLVE_UPD', 'ITER_ASSM', label='ITER = ITER + 1', fontname='sans-serif', fontsize='9')
        c3.edge('CONV_CHK', 'ACCUM_PL', label='Yes (Converged)', fontname='sans-serif', fontsize='9')
        c3.edge('ACCUM_PL', 'MORE_INC')
        c3.edge('MORE_INC', 'INC_START', label='Yes (FMULT = FMULT · FINCR)', fontname='sans-serif', fontsize='9')
        c3.edge('MORE_INC', 'RECONSTR', label='No (All Increments Done)', fontname='sans-serif', fontsize='9')

    # Global connections between phases
    dot.edge('START', 'INT_CELL')
    dot.edge('EVAL_KE', 'MAX_YIELD')
    dot.edge('DEC_ELAST', 'END_NL', label='Yes (Skip Plasticity)', fontname='sans-serif', fontsize='9')
    dot.edge('FORM_MAT', 'INC_START')
    dot.edge('RECONSTR', 'END_NL')

    # Legend
    with dot.subgraph(name='cluster_nl_legend') as leg:
        leg.attr(label='Legend', fontname='sans-serif', fontsize='10.5', color='black', style='solid', bgcolor='#ffffff')
        leg.node('l1', '2D Cell Integration', **cell_attr)
        leg.node('l2', 'Increment / State Update', **plast_attr)
        leg.node('l3', 'Inner Multiplier Iteration', **iter_attr)
        leg.node('l4', 'Decision / Convergence', **dec_attr)
        leg.edge('l1', 'l2', style='invis')
        leg.edge('l2', 'l3', style='invis')
        leg.edge('l3', 'l4', style='invis')

    dot.render('diagrama_nl_graphviz', format='png', cleanup=True)
    print("O diagrama detalhado da parte não-linear foi exportado com sucesso como 'diagrama_nl_graphviz.png'!")


if __name__ == '__main__':
    build_global_flowchart()
    build_nonlinear_flowchart()
