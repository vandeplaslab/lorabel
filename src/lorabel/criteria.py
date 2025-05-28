import math
from typing import NoReturn, Union

from ..math import *


class CRITERIA:
    """The criteria class provides the processing and checking of criteria for decompositions.

    Parameters
    ----------
    criteria : list
        List containing criteria, each item (a.k.a. criteria) can be a list or a simple string, e.g.
            [['$a$', '<', '0.05']] or [['$a$', '<', '0.05'], ['linalg.norm($a$,2)+linalg.norm($b$, 2)', '<', '0.05']].
        variables need to be put between dollar signs, i.e. $.$.

    condition : {'any', 'all'}, default='all'
        Condition for the list of criteria, for 'any' check returns `True`
        if one or more critera are met. If 'all', check returns `True` if all critera are met.

    dict_x : dict, default={}
        Dictionary linking $.$ variables to variables in the x dictionary providing the actual variables.

    state : list
        Current state of the criteria, returns list with left and right parts (0 and 2 indices) of the condition

    Attributes
    ----------
    check : Check if criteria are fulfilled

    _evaluate : Evaluate a single criterium

    _stringify : Go back to string version and join list to single string

    _format_list : Format the list of criteria

    _create_dictx : Creates a dictionary for x to replace $.$ variables in the criteria

    _replace_dict : Replaces $.$ variables with the corresponding variables as defined by _create_dictx

    """

    def __init__(
        self,
        criteria: list,
        condition: str = "all",
    ) -> NoReturn:
        self.criteria = [criteria] if isinstance(criteria[0], str) is True else criteria
        self.condition = eval(condition)
        self.dict_x = {}
        self.state = []

        pass

    def check(self, x: dict) -> bool:
        """Check if criteria are fulfilled or returns True if no criteria are set

        Parameters
        ----------
        x : dict
            Dictionary linking variables to memory ids

        Returns
        -------
            bool : True if conditions are fulfilled, false if they are not
        """
        if self.criteria == [[""]]:
            return True
        else:
            self.dict_x = self._create_dictx(x)
            self.criteria_list = self._format_list()
            self.state = [
                self._evaluate(criteria, x) for criteria in self.criteria_list
            ]
            list_of_strings = [
                self._stringify(single_state) for single_state in self.state
            ]
            check = eval(
                "condition(criteria_list)",
                {
                    "condition": self.condition,
                    "criteria_list": [
                        eval(string_evaluate) for string_evaluate in list_of_strings
                    ],
                },
            )

            return check

    def _stringify(self, a: list) -> str:
        """Go back to string version and join list to single string

        Parameters
        ----------
        a : list
            List of strings containing a single criterium that has been evaluated

        Returns
        -------
            str : String containing a single condition that has been evaluated
        """
        a = [str(a_part) for a_part in a]

        return "".join(a)

    def _evaluate(self, a: list[str, str, str], x: dict) -> list:
        """Evaluate a single criterium

        Parameters
        ----------
        a : list[str, str, str]
            List of strings containing a single criterium

        x : dict
            Dictionary containing values

        Returns
        -------
            list : List with eval() values.
        """
        a[0] = eval(a[0], None, {"x": x, "linalg": linalg, "math": math})
        a[2] = eval(a[2], None, {"x": x, "linalg": linalg, "math": math})

        return a

    def _format_list(
        self,
    ) -> list[list]:
        """Format the list of criteria, replacing $xi$ values by dict_x

        Returns
        -------
            list[list] : Generate list with criteria that consists of list of strings
        """
        criteria_list = [self._replace_dict(criteria) for criteria in self.criteria]

        return criteria_list

    def _create_dictx(
        self,
        x: dict,
    ) -> dict:
        """Creates a dictionary for x to replace $.$ variables in the criteria

        Parameters
        ----------
        x : dict
            Dictionary linking variables to memory ids

        Returns
        -------
            dict : dictionary linking $.$ to the actual values of the x dict
        """
        dict_x = {"$" + a + "$": 'x["' + a + '"]' for a in x.keys()}

        return dict_x

    def _replace_dict(
        self,
        criteria: list[str, str, str],
    ) -> list[str, str, str]:
        """Replaces $.$ variables with the corresponding variables as defined by _create_dictx

        Parameters
        ----------
        criteria : list[str, str, str]
            List containing a single criteria as three strings

        Returns
        -------
            list[str, str, str] : List of strings where $.$ variables are replaced by proper variables in the x dictionary

        """
        for key in self.dict_x.keys():
            criteria = [subpart.replace(key, self.dict_x[key]) for subpart in criteria]

        return criteria
