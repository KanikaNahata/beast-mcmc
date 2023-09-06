/*
 * GeneralizedAdditiveGaussianProcessModelParser.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inferencexml.distribution;

import dr.inference.distribution.*;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.inferencexml.distribution.GeneralizedLinearModelParser.*;
import static dr.inferencexml.glm.GeneralizedLinearModelParser.DEPENDENT_VARIABLES;

/**
 * @author Filippo Monti
 * @author Marc Suchard
 */
public class GeneralizedAdditiveGaussianProcessModelParser extends AbstractXMLObjectParser {

    public static final String GAM_GP_LIKELIHOOD = "gamGpModel";

    public String getParserName() {
        return GAM_GP_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(DEPENDENT_VARIABLES);
        Parameter dependentParam = null;
        if (cxo != null) {
            dependentParam = (Parameter) cxo.getChild(Parameter.class);
        }

//        String family = xo.getStringAttribute(FAMILY);
        LogGaussianProcessModel gp = new LogGaussianProcessModel(dependentParam);

        addIndependentParameters(xo, gp, dependentParam);

        return gp;
    }

//    public void addRandomEffects(XMLObject xo, GeneralizedLinearModel glm,
//                                 Parameter dependentParam) throws XMLParseException {
//        int totalCount = xo.getChildCount();
//
//        for (int i = 0; i < totalCount; i++) {
//            if (xo.getChildName(i).compareTo(RANDOM_EFFECTS) == 0) {
//                XMLObject cxo = (XMLObject) xo.getChild(i);
//                Parameter randomEffect = (Parameter) cxo.getChild(Parameter.class);
//                checkRandomEffectsDimensions(randomEffect, dependentParam);
//                glm.addRandomEffectsParameter(randomEffect);
//            }
//        }
//    }

    // TODO remove code duplication
    public void addIndependentParameters(XMLObject xo, GeneralizedLinearModel glm,
                                         Parameter dependentParam) throws XMLParseException {
        int totalCount = xo.getChildCount();

        for (int i = 0; i < totalCount; i++) {
            if (xo.getChildName(i).compareTo(INDEPENDENT_VARIABLES) == 0) {
                XMLObject cxo = (XMLObject) xo.getChild(i);
                Parameter independentParam = (Parameter) cxo.getChild(Parameter.class);
                DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
                checkDimensions(independentParam, dependentParam, designMatrix);
                cxo = cxo.getChild(INDICATOR);
                Parameter indicator = null;
                if (cxo != null) {
                    indicator = (Parameter) cxo.getChild(Parameter.class);
                    if (indicator.getDimension() <= 1) {
                        // if a dimension hasn't been set, then set it automatically
                        indicator.setDimension(independentParam.getDimension());
                    }
                    if (indicator.getDimension() != independentParam.getDimension())
                        throw new XMLParseException("dim(" + independentParam.getId() + ") != dim(" + indicator.getId() + ")");
                }

//                if (checkFullRankOfMatrix) {
//                    checkFullRank(designMatrix);
//                }

                glm.addIndependentParameter(independentParam, designMatrix, indicator);
            }
        }
    }

//    private boolean checkFullRankOfMatrix = false;

//    private void checkFullRank(DesignMatrix designMatrix) throws XMLParseException {
//        int fullRank = designMatrix.getColumnDimension();
////        System.err.println("designMatrix getColumnDimension = "+fullRank);
//        SingularValueDecomposition svd = new SingularValueDecomposition(
//                new DenseDoubleMatrix2D(designMatrix.getParameterAsMatrix()));
//        int realRank = svd.rank();
//        if (realRank != fullRank) {
//            throw new XMLParseException(
//                "rank(" + designMatrix.getId() + ") = " + realRank +
//                        ".\nMatrix is not of full rank as colDim(" + designMatrix.getId() + ") = " + fullRank
//            );
//        }
//    }

//    private void checkRandomEffectsDimensions(Parameter randomEffect, Parameter dependentParam)
//            throws XMLParseException {
//        if (dependentParam != null) {
//            if (randomEffect.getDimension() <= 1) {
//                // if a dimension hasn't been set, then set it automatically
//                randomEffect.setDimension(dependentParam.getDimension());
//            }
//            if (randomEffect.getDimension() != dependentParam.getDimension()) {
//                throw new XMLParseException(
//                        "dim(" + dependentParam.getId() + ") != dim(" + randomEffect.getId() + ")"
//                );
//            }
//        }
//    }

    // TODO remove code duplication
    private void checkDimensions(Parameter independentParam, Parameter dependentParam, DesignMatrix designMatrix)
            throws XMLParseException {
        if (dependentParam != null) {
            if (dependentParam.getDimension() <= 1) {
                // if a dimension hasn't been set, then set it automatically
                dependentParam.setDimension(designMatrix.getRowDimension());
            }
            if ((dependentParam.getDimension() != designMatrix.getRowDimension()) ||
                    (independentParam.getDimension() != designMatrix.getColumnDimension()))
                throw new XMLParseException(
                        "dim(" + dependentParam.getId() + ") != dim(" + designMatrix.getId() + " %*% " + independentParam.getId() + ")"
                );
        } else {
            if (independentParam.getDimension() <= 1) {
                // if a dimension hasn't been set, then set it automatically
                independentParam.setDimension(designMatrix.getColumnDimension());
            }
            if (independentParam.getDimension() != designMatrix.getColumnDimension()) {
                throw new XMLParseException(
                        "dim(" + independentParam.getId() + ") is incompatible with dim (" + designMatrix.getId() + ")"
                );
            }
//            System.err.println(independentParam.getId()+" and "+designMatrix.getId());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            AttributeRule.newStringRule(FAMILY),
//            AttributeRule.newBooleanRule(CHECK_IDENTIFIABILITY, true),
//            AttributeRule.newBooleanRule(CHECK_FULL_RANK, true),
            new ElementRule(DEPENDENT_VARIABLES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(INDEPENDENT_VARIABLES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, true),
                            new ElementRule(DesignMatrix.class),
                            new ElementRule(INDICATOR,
                                    new XMLSyntaxRule[]{
                                            new ElementRule(Parameter.class)
                                    }, true),
                    }, 0, 10),
//            new ElementRule(RANDOM_EFFECTS,
//                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, 0, 3),
    };

    public String getParserDescription() {
        return "Calculates the generalized linear model likelihood of the dependent parameters given one or more blocks of independent parameters and their design matrix.";
    }

    public Class getReturnType() {
        return LogGaussianProcessModel.class;
    }
}
