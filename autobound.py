# ------------------------------------------------------------
# Bounds for ATE, CDE, and NDE using Autobounds
# ------------------------------------------------------------
# Variables:
#   z : instrument / protected attribute
#   x : mediator / treatment
#   y : outcome
#
# This script:
# 1. loads observed joint probabilities for (z, x, y),
# 2. defines two DAGs:
#       - standard IV model
#       - IV model with exclusion restriction violation (z -> y)
# 3. computes ATE bounds under both models,
# 4. computes CDE(0), CDE(1), NDE(0), NDE(1) bounds
#    under the exclusion-restriction-violated model.
# ------------------------------------------------------------

from autobounds.causalProblem import causalProblem
from autobounds.DAG import DAG
import pandas as pd


# ------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------

def make_data() -> pd.DataFrame:
    """Create the observed joint distribution P(z, x, y)."""
    data = pd.DataFrame([
        {"z": 0, "x": 1, "y": 1, "n": 527},
        {"z": 0, "x": 1, "y": 0, "n": 527},
        {"z": 0, "x": 0, "y": 1, "n": 55},
        {"z": 0, "x": 0, "y": 0, "n": 303},
        {"z": 1, "x": 1, "y": 1, "n": 880},
        {"z": 1, "x": 1, "y": 0, "n": 195},
        {"z": 1, "x": 0, "y": 1, "n": 5},
        {"z": 1, "x": 0, "y": 0, "n": 59},
    ])

    data["prob"] = data["n"] / data["n"].sum()
    return data.drop(columns=["n"])


# ------------------------------------------------------------
# 2. DAG definitions
# ------------------------------------------------------------

def make_standard_iv_dag() -> DAG:
    """
    Standard IV DAG:
        z -> x -> y
        u -> x
        u -> y
    """
    return DAG(
        edges="z -> x, x -> y, u -> x, u -> y",
        unob="u"
    )


def make_er_violated_dag() -> DAG:
    """
    IV DAG with exclusion restriction violation:
        z -> x -> y
        z -> y
        u -> x
        u -> y
    """
    return DAG(
        edges="z -> x, x -> y, z -> y, u -> x, u -> y",
        unob="u"
    )


# ------------------------------------------------------------
# 3. Generic solver helpers
# ------------------------------------------------------------

def make_problem(dag: DAG, data: pd.DataFrame) -> causalProblem:
    """Initialize a causal problem and load observed data."""
    prob = causalProblem(dag)
    prob.load_data(data)
    return prob


def solve_ate(dag: DAG, data: pd.DataFrame):
    """Solve ATE bounds for a given DAG."""
    prob = make_problem(dag, data)
    prob.set_ate(ind="x", dep="y")
    return prob.solve()


def solve_custom_estimand(dag: DAG, data: pd.DataFrame, estimand_builder):
    """
    Solve bounds for a custom estimand.

    Parameters
    ----------
    dag : DAG
        Causal graph.
    data : pd.DataFrame
        Observed probability table.
    estimand_builder : callable
        Function taking `prob` and returning an estimand expression.
    """
    prob = make_problem(dag, data)
    prob.set_estimand(estimand_builder(prob))
    return prob.solve()


# ------------------------------------------------------------
# 4. Estimand builders
# ------------------------------------------------------------

def cde0(prob):
    """
    CDE(0) = E[Y(x=0, z=1) - Y(x=0, z=0)]
    """
    return prob.p("y(x=0,z=1)=1") - prob.p("y(x=0,z=0)=1")


def cde1(prob):
    """
    CDE(1) = E[Y(x=1, z=1) - Y(x=1, z=0)]
    """
    return prob.p("y(x=1,z=1)=1") - prob.p("y(x=1,z=0)=1")


def nde0(prob):
    """
    NDE(0) = E[Y(X(z=0), z=1) - Y(X(z=0), z=0)]

    Expanded as:
      P(Y(x=0,z=1)=1, X(z=0)=0) - P(Y(x=0,z=0)=1, X(z=0)=0)
    + P(Y(x=1,z=1)=1, X(z=0)=1) - P(Y(x=1,z=0)=1, X(z=0)=1)
    """
    term0 = (
        prob.p("y(x=0,z=1)=1&x(z=0)=0")
        - prob.p("y(x=0,z=0)=1&x(z=0)=0")
    )
    term1 = (
        prob.p("y(x=1,z=1)=1&x(z=0)=1")
        - prob.p("y(x=1,z=0)=1&x(z=0)=1")
    )
    return term0 + term1


def nde1(prob):
    """
    NDE(1) = E[Y(X(z=1), z=1) - Y(X(z=1), z=0)]

    Expanded as:
      P(Y(x=0,z=1)=1, X(z=1)=0) - P(Y(x=0,z=0)=1, X(z=1)=0)
    + P(Y(x=1,z=1)=1, X(z=1)=1) - P(Y(x=1,z=0)=1, X(z=1)=1)
    """
    term0 = (
        prob.p("y(x=0,z=1)=1&x(z=1)=0")
        - prob.p("y(x=0,z=0)=1&x(z=1)=0")
    )
    term1 = (
        prob.p("y(x=1,z=1)=1&x(z=1)=1")
        - prob.p("y(x=1,z=0)=1&x(z=1)=1")
    )
    return term0 + term1


# ------------------------------------------------------------
# 5. Main execution
# ------------------------------------------------------------

def main():
    data = make_data()

    dag_standard = make_standard_iv_dag()
    dag_er_violated = make_er_violated_dag()

    print("=== Standard IV model ===")
    result_standard = solve_ate(dag_standard, data)
    print("ATE bounds:", result_standard)

    print("\n=== Exclusion restriction violated model ===")
    result_violated = solve_ate(dag_er_violated, data)
    print("ATE bounds:", result_violated)

    print("\n=== ER-violated model: direct-effect bounds ===")
    print("CDE0 bounds:", solve_custom_estimand(dag_er_violated, data, cde0))
    print("CDE1 bounds:", solve_custom_estimand(dag_er_violated, data, cde1))
    print("NDE0 bounds:", solve_custom_estimand(dag_er_violated, data, nde0))
    print("NDE1 bounds:", solve_custom_estimand(dag_er_violated, data, nde1))


if __name__ == "__main__":
    main()