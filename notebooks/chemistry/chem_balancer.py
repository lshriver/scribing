# chem_balancer.py

import numpy as np
from fractions import Fraction
from sympy import Matrix, lcm
import re

def parse_compound(compound):
    """
    Parses a chemical compound with support for nested parentheses.
    E.g., "Al2(SO4)3" -> {'Al': 2, 'S': 3, 'O': 12}
    """
    pattern = r'([A-Z][a-z]?|$$|$$|\d+)'
    tokens = re.findall(pattern, compound)
    def parse_tokens(tokens):
        stack = [{}]
        i = 0
        while i < len(tokens):
            token = tokens[i]
            if token == '(':
                stack.append({})
                i += 1
            elif token == ')':
                group = stack.pop()
                i += 1
                multiplier = 1
                if i < len(tokens) and tokens[i].isdigit():
                    multiplier = int(tokens[i])
                    i += 1
                for element, count in group.items():
                    stack[-1][element] = stack[-1].get(element, 0) + count * multiplier
            elif re.match(r'[A-Z][a-z]?$', token):
                element = token
                count = 1
                i += 1
                if i < len(tokens) and tokens[i].isdigit():
                    count = int(tokens[i])
                    i += 1
                stack[-1][element] = stack[-1].get(element, 0) + count
            else:
                # It's a digit, already handled
                i += 1
        return stack[0]
    return parse_tokens(tokens)

def build_matrix(equation):
    """
    Builds a stoichiometry matrix from the chemical equation.
    """
    reactants_str, products_str = equation.replace(' ', '').split('=')
    reactants = reactants_str.split('+')
    products = products_str.split('+')

    all_compounds = reactants + products
    element_list = []
    compound_elements = []

    # Parse compounds and collect elements
    for compound in all_compounds:
        elements = parse_compound(compound)
        compound_elements.append(elements)
        for el in elements:
            if el not in element_list:
                element_list.append(el)

    # Build the matrix
    num_elements = len(element_list)
    num_compounds = len(all_compounds)

    matrix = np.zeros((num_elements, num_compounds), dtype=int)

    for i, element in enumerate(element_list):
        for j, compound_dict in enumerate(compound_elements):
            count = compound_dict.get(element, 0)
            if j < len(reactants):
                matrix[i][j] = count
            else:
                # Products have negative counts
                matrix[i][j] = -count

    return matrix, all_compounds, element_list

def balance_equation(equation):
    """
    Balances the chemical equation and returns it as a string.
    """
    matrix, compounds, elements = build_matrix(equation)
    num_elements, num_compounds = matrix.shape

    # Remove one equation (make the system underdetermined)
    matrix_reduced = matrix[:-1, :]  # Remove the last row

    # Convert to sympy Matrix for exact calculations
    sympy_matrix = Matrix(matrix_reduced)
    null_space = sympy_matrix.nullspace()

    if not null_space:
        raise ValueError("No solution found. The equation may already be balanced or invalid.")

    basis = null_space[0]
    lcm_denoms = lcm([term.q for term in basis])
    coeffs = [term * lcm_denoms for term in basis]
    coeffs = [int(c) for c in coeffs]
    coeffs = normalize_coefficients(coeffs)

    # Build the balanced equation string
    num_reactants = len(equation.replace(' ', '').split('=')[0].split('+'))
    balanced_reactants = []
    for i in range(num_reactants):
        coeff = coeffs[i]
        compound = compounds[i]
        if coeff != 1:
            balanced_reactants.append(f'{coeff}{compound}')
        else:
            balanced_reactants.append(compound)
    balanced_products = []
    for i in range(num_reactants, len(compounds)):
        coeff = coeffs[i]
        compound = compounds[i]
        if coeff != 1:
            balanced_products.append(f'{coeff}{compound}')
        else:
            balanced_products.append(compound)
    balanced_eq = ' + '.join(balanced_reactants) + ' = ' + ' + '.join(balanced_products)
    return balanced_eq

def normalize_coefficients(coeffs):
    """
    Normalizes coefficients to the smallest integer values.
    """
    gcd_coeffs = np.gcd.reduce(coeffs)
    coeffs = [c // gcd_coeffs for c in coeffs]
    # Ensure all coefficients are positive
    if any(c < 0 for c in coeffs):
        coeffs = [-c for c in coeffs]
    return coeffs

# If the module is run as a script, allow input from the user
if __name__ == "__main__":
    equation = input("Enter the unbalanced chemical equation:\n")
    try:
        balanced_equation = balance_equation(equation)
        print('Balanced Equation:')
        print(balanced_equation)
    except Exception as e:
        print(f'Error: {e}')