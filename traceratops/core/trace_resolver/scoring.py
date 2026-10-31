"""Orthogonal history selection and empirical distance costs."""

import math

import numpy as np


class PolymerScorer:
    def __init__(
        self,
        distance_model,
        history_mode="multi",
        distance_score="residual",
        history_length=3,
    ):
        if history_mode not in {"nearest", "multi"}:
            raise ValueError("history_mode must be 'nearest' or 'multi'")
        if distance_score not in {"residual", "likelihood"}:
            raise ValueError("distance_score must be 'residual' or 'likelihood'")
        if history_length < 1:
            raise ValueError("history_length must be positive")
        self.model = distance_model
        self.history_mode = history_mode
        self.distance_score = distance_score
        self.history_length = history_length

    def transition_cost(self, coordinate, position, history):
        if not history:
            return 0.0
        previous = (
            history[-1:]
            if self.history_mode == "nearest"
            else history[-self.history_length :]
        )
        costs = []
        for old_position, old_coordinate in previous:
            distance = float(np.linalg.norm(coordinate - old_coordinate))
            stats = self.model.for_separation(position - old_position)
            residual = (distance - stats.mean) / stats.standard_deviation
            if self.distance_score == "residual":
                costs.append(residual * residual)
            else:
                costs.append(
                    0.5 * residual * residual
                    + math.log(stats.standard_deviation)
                    + 0.5 * math.log(2 * math.pi)
                )
        return float(np.mean(costs))
